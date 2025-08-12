# === Load Libraries and Sources ===
suppressPackageStartupMessages({
library(optparse)
library(Biostrings)
library(ggplot2)
library(patchwork)
})

source("~/Documents/basic_gyrae/scripts/handle_EVO2_output.r")
source("~/Documents/basic_gyrae/scripts/evo2_analysis_functions.R")

# === Command-Line Interface ===
option_list <- list(
  make_option(c("-f", "--fasta"), type = "character", help = "FASTA file path"),
  make_option(c("-v", "--variants"), type = "character", help = "Comma-separated variant names (e.g., 83_S1,83_L,...)"),
  make_option(c("--highlight"), type = "character", default = "247,248,249", help = "Comma-separated highlight positions"),
  make_option(c("-i", "--input_logits_dir"), type = "character", help = "path to input logits directory"),
  make_option(c("-o", "--output_dir"), type = "character", help = "path to output directory"),
  make_option(c("-r", "--reverse_comp"), type = "logical", action = "store_true", default = FALSE,
              help = "Use reverse complement sequence"),
  make_option(c("-s", "--only_single"), type = "logical", action = "store_true", default = FALSE,
              help = "Only save single plots"),
  make_option(c("--index_start"), type = "integer", default = 200,
              help = "Start position for plotting range"),
  make_option(c("--index_end"), type = "integer", default = 300,
              help = "End position for plotting range")
)

opt <- parse_args(OptionParser(option_list = option_list))
fasta_file <- opt$fasta
variant_names <- unlist(strsplit(opt$variants, ","))

highlight <- as.integer(unlist(strsplit(opt$highlight, ",")))

reverse_comp <- opt$reverse_comp
only_single <- opt$only_single

input_dir  <- opt$input_logits_dir
output_dir <- opt$output_dir

# Create index range for plotting
index_rows <- opt$index_start:opt$index_end
cat("\nPlotting range:", min(index_rows), "to", max(index_rows), "\n")

# === Load Sequences ===
all_variants <- readDNAStringSet(fasta_file)

# === Initialize DataFrames ===
variant_dfs <- list()
for (name in variant_names) {
  npy_file <- paste0(input_dir, "/input_", name, "_logits.npy")
  variant_dfs[[name]] <- initialize_df(npy_file, name, all_variants)
}

# === Plot Entropy & Log-Likelihood ===
entropy_plots <- list()
ll_plots <- list()
stacked_single <- list()
accuracy_results <- list()
for (name in names(variant_dfs)) {
  df <- variant_dfs[[name]]
  entropy_plots[[name]] <- plot_entropy(name, df, 
                                      index_rows = index_rows,
                                      highlighted = highlight, 
                                      reverse_comp = reverse_comp)
  ll_plots[[name]] <- plot_log_likelihood(name, df, 
                                         index_rows = index_rows,
                                         highlighted = highlight, 
                                         reverse_comp = reverse_comp)
  # Stack entropy and log-likelihood vertically for concise single plot
  # Strip x-axis on top plot to remove fluff when stacked
  stacked_single[[name]] <- ((entropy_plots[[name]] + strip_x()) / ll_plots[[name]]) +
    plot_layout(heights = c(1, 1), guides = "collect") &
    theme(legend.position = "right")
  accuracy_results[[name]] <- probability_correct_base(name, df, all_variants)
}

cat("\n✅ Prediction Accuracy:\n")
for (result in accuracy_results) {
  cat(result, "\n")
}

# === Plot Differences vs First Variant (as reference) ===

# treating first variant name as reference for compare
ref_name <- variant_names[1]
ref_df <- variant_dfs[[ref_name]]

diff_plots <- list()
entropy_diffs <- list()
loglike_diffs <- list()

# variant_names[-1] skips first element (reference gene)
for (name in variant_names[-1]) {
  # store ref df
  target_df <- variant_dfs[[name]]
  entropy_diffs[[name]] <- target_df$entropy - ref_df$entropy
  loglike_diffs[[name]] <- target_df$log_likelihood - ref_df$log_likelihood
  
  diff_plots[[paste0("entropy_", name, "-", ref_name)]] <- plot_metric_diff(ref_name, name, ref_df, target_df, 
                                    metric = "entropy", 
                                    index_rows = index_rows,
                                    highlight = highlight, 
                                    reverse_comp = reverse_comp)
  diff_plots[[paste0("ll_", name, "-", ref_name)]] <- plot_metric_diff(ref_name, name, ref_df, target_df, 
                                    metric = "log_likelihood", 
                                    index_rows = index_rows,
                                    highlight = highlight, 
                                    reverse_comp = reverse_comp)
}

# Build stacked difference figures: all entropy diffs stacked; all ll diffs stacked
entropy_diff_keys <- grep("^entropy_", names(diff_plots), value = TRUE)
ll_diff_keys <- grep("^ll_", names(diff_plots), value = TRUE)

stacked_entropy_diffs <- NULL
stacked_ll_diffs <- NULL
if (length(entropy_diff_keys) > 0) {
  # Apply strip_x to all but the bottom plot
  entropy_stack_list <- diff_plots[entropy_diff_keys]
  if (length(entropy_stack_list) > 1) {
    for (i in seq_len(length(entropy_stack_list) - 1)) {
      entropy_stack_list[[i]] <- entropy_stack_list[[i]] + strip_x()
    }
  }
  stacked_entropy_diffs <- wrap_plots(entropy_stack_list, ncol = 1, guides = "collect") &
    theme(legend.position = "right")
}
if (length(ll_diff_keys) > 0) {
  ll_stack_list <- diff_plots[ll_diff_keys]
  if (length(ll_stack_list) > 1) {
    for (i in seq_len(length(ll_stack_list) - 1)) {
      ll_stack_list[[i]] <- ll_stack_list[[i]] + strip_x()
    }
  }
  stacked_ll_diffs <- wrap_plots(ll_stack_list, ncol = 1, guides = "collect") &
    theme(legend.position = "right")
}

# === Log-Likelihood Matrix & Heatmap ===
total_ll <- sapply(variant_dfs, function(df) sum(df$log_likelihood))
diff_matrix <- outer(total_ll, total_ll, "-")
df_plot <- as.data.frame(as.table(diff_matrix))

heatmap_plot <- ggplot(df_plot, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", Freq))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Δ log-L") +
  labs(title = "Log-Likelihood Differences Between Variants (col-row)", x = "", y = "") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

# === Save All Plots ===
save_pdfs <- function(output_dir, only_single = FALSE) {
  # Save single plots
  for (name in names(stacked_single)) {
    # Height scaled for two rows
    save_as_pdf(stacked_single[[name]], paste0("", name, ".pdf"), out_dir = paste0(output_dir,"/single"), height = 10)
  }
  # Save heatmap
    save_as_pdf(heatmap_plot, "tot_ll.pdf", out_dir = paste0(output_dir, "/tot_ll"))
  
  # Always save stacked difference figures (entropy and log-likelihood)
  if (!is.null(stacked_entropy_diffs)) {
    # scale height with number of panels (approx 3 inches per panel)
    save_as_pdf(stacked_entropy_diffs, "entropy_diffs_stacked.pdf", out_dir = paste0(output_dir, "/compare"), height = max(6, 3 * length(entropy_diff_keys)))
  }
  if (!is.null(stacked_ll_diffs)) {
    save_as_pdf(stacked_ll_diffs, "ll_diffs_stacked.pdf", out_dir = paste0(output_dir, "/compare"), height = max(6, 3 * length(ll_diff_keys)))
  }
}

# Save plots with the only_single flag
save_pdfs(output_dir, only_single = only_single)

cat("\n✅ All plots generated.")