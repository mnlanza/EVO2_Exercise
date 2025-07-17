# === Load Libraries and Sources ===
suppressPackageStartupMessages({
library(optparse)
library(Biostrings)
library(ggplot2)
})

source("handle_EVO2_output.r")
source("evo2_analysis_functions.R")

# === Command-Line Interface ===
option_list <- list(
  make_option(c("-f", "--fasta"), type = "character", help = "FASTA file path"),
  make_option(c("-v", "--variants"), type = "character", help = "Comma-separated variant names (e.g., 83_S1,83_L,...)"),
  make_option(c("--highlight"), type = "character", default = "247,248,249", help = "Comma-separated highlight positions"),
  make_option(c("-i", "--input_dir"), type = "character", help = "path to input directory"),
  make_option(c("-o", "--output_dir"), type = "character", help = "path to output directory"),
  make_option(c("-r", "--reverse_comp"), type = "logical", action = "store_true", default = FALSE,
              help = "Use reverse complement sequence")
)

opt <- parse_args(OptionParser(option_list = option_list))
fasta_file <- opt$fasta
variant_names <- unlist(strsplit(opt$variants, ","))

highlight <- as.integer(unlist(strsplit(opt$highlight, ",")))

reverse_comp <- opt$reverse_comp

input_dir  <- opt$input_dir
output_dir <- opt$output_dir


###TESTING####
# setwd("~/Desktop/Relman_lab/EVO2_Exercise")
# file.exists("~/Desktop/Relman_lab/All_variant_likelihoods/reverse_complements_83.fasta")
# fasta_file <- "~/Desktop/Relman_lab/All_variant_likelihoods/reverse_complements_83.fasta"
# highlight <- 247:249
# variant_names <- c("83_AGC_rc", "83_TGA_rc", "83_AGT_rc", "83_ATC_rc")
# reverse_comp <- TRUE
# input_dir <- "~/Desktop/Relman_lab/All_variant_likelihoods/raw_data"
# file.exists("~/Desktop/Relman_lab/All_variant_likelihoods/raw_data/input_83_AGC_rc_logits.npy")

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
accuracy_results <- list()
for (name in names(variant_dfs)) {
  df <- variant_dfs[[name]]
  entropy_plots[[name]] <- plot_entropy(name, df, highlighted = highlight, reverse_comp = reverse_comp)
  ll_plots[[name]] <- plot_log_likelihood(name, df, highlighted = highlight, reverse_comp = reverse_comp)
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
                                    metric = "entropy", highlight = highlight, reverse_comp = reverse_comp)
  diff_plots[[paste0("ll_", name, "-", ref_name)]] <- plot_metric_diff(ref_name, name, ref_df, target_df, 
                                    metric = "log_likelihood", highlight = highlight, reverse_comp = reverse_comp)
}

# === Log-Likelihood Matrix & Heatmap ===
total_ll <- sapply(variant_dfs, function(df) sum(df$log_likelihood))
diff_matrix <- outer(total_ll, total_ll, "-")
df_plot <- as.data.frame(as.table(diff_matrix))

heatmap_plot <- ggplot(df_plot, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", Freq))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Δ log-L") +
  labs(title = "Log-Likelihood Differences Between Variants", x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# === Save All Plots ===
save_pdfs <- function(output_dir) {
  # Save single plots
  for (name in names(entropy_plots)) {
    save_as_pdf(entropy_plots[[name]], paste0("entropy_", name, ".pdf"), out_dir = paste0(output_dir,"/single"))
    save_as_pdf(ll_plots[[name]], paste0("ll_", name, ".pdf"), out_dir= paste0(output_dir,"/single"))
  }
  
  # Save diff plots
  for (name in names(diff_plots)) {
    save_as_pdf(diff_plots[[name]], paste0(name, ".pdf"), out_dir = paste0(output_dir,"/compare"))
  }
  
  # Save heatmap
  save_as_pdf(heatmap_plot, "tot_ll.pdf", out_dir = paste0(output_dir, "/tot_ll"))
}

# Uncomment this to save plots
save_pdfs(output_dir)

cat("\n✅ All plots generated.")