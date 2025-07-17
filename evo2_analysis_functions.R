#source("handle_EVO2_output.r")
suppressPackageStartupMessages({
library(Biostrings)
library(ggplot2)
})
#' Initialize dataframe with EVO2 probabilities, entropy, and log-likelihood

#' @param EVO2_npy_file Path to the numpy file containing logits
#' @param sequence_name Name of the sequence to analyze
#' @param all_variants DNAStringSet object containing the variants
#' @return Dataframe with probabilities, entropy, and log-likelihood

initialize_df <- function(EVO2_npy_file, sequence_name, all_variants) { 
  # Load and normalize logits into nucleotide probabilities
  prob_matrix <- get_logits(EVO2_npy_file) 

  # Calculate entropy and add it as a column
  prob_matrix$entropy <- apply(prob_matrix[,1:4], 1, 
    function(position_prob) -sum(position_prob * log2(position_prob)))
  # Get sequence and convert to nucleotide list
  seqs <- as.character(all_variants)
  desired_seq <- seqs[sequence_name] 
  seq_vec <- strsplit(desired_seq, split = "")[[1]] 
  
  # Calculate log-likelihood
  log_likelihood <- mapply(function(base, i) {
      base_index <- match(base, c("A", "C", "G", "T"))
      log2(prob_matrix[i, base_index])}, seq_vec, seq_along(seq_vec))
  prob_matrix$log_likelihood <- log_likelihood
  return(prob_matrix)
}

#' Build dataframe for plotting
#' 
#' @param gene_df Source dataframe with metrics
#' @param metric Name of metric column to plot
#' @param index_rows Range of positions to plot
#' @param highlighted Positions to highlight
#' @return Dataframe ready for plotting
build_plot_df <- function(gene_df, metric, index_rows = 200:300, highlighted = 247:249) {
  df <- data.frame(
    pos     = index_rows,
    value   = gene_df[[metric]][index_rows],
    is_high = index_rows %in% highlighted,
    codon_i = ((index_rows - 1) %% 3) + 1
  )
  df
}

#' Plot entropy across positions
#' 
#' @param gene_name Name of the gene variant
#' @param gene_df Dataframe containing entropy data
#' @param index_rows Range of positions to plot
#' @param highlighted Positions to highlight
#' @return ggplot object
plot_entropy <- function(gene_name, gene_df, index_rows = 200:300, highlighted = 247:249, reverse_comp = FALSE) {
  if (reverse_comp){
    highlighted <- (nrow(gene_df) - highlighted + 1)
    index_rows <- sort(nrow(gene_df) - index_rows + 1)
  }
  plot_df <- build_plot_df(gene_df, "entropy", index_rows, highlighted)
  ggplot(plot_df, aes(x = pos, y = value)) + 
    geom_point(aes(color = factor(codon_i), size = is_high),
               show.legend = c(color = TRUE, size = TRUE)) + 
    scale_color_manual(name = "Codon position",
                      values = c("1" = "red", "2" = "orange", "3" = "green")) +
    scale_size_manual(name = "Key positions",
                     values = c(`FALSE` = 1.5, `TRUE` = 2.7)) +
    labs(x = "Nucleotide position", y = "Entropy",
         title = sprintf("Entropy across positions %d-%d (%s)", 
                        min(index_rows), max(index_rows), gene_name)) +
    theme_minimal() + theme(legend.position = "right", legend.box = "vertical")
}

#' Plot log-likelihood across positions
#' 
#' @param gene_name Name of the gene variant
#' @param gene_df Dataframe containing log-likelihood data
#' @param index_rows Range of positions to plot
#' @param highlighted Positions to highlight
#' @return ggplot object
plot_log_likelihood <- function(gene_name, gene_df, index_rows = 200:300, highlighted = 247:249, reverse_comp = FALSE) {
  if (reverse_comp){
    highlighted <- (nrow(gene_df) - highlighted + 1)
    index_rows <- sort(nrow(gene_df) - index_rows + 1)
  }
  plot_df <- build_plot_df(gene_df, "log_likelihood", index_rows, highlighted)
  ggplot(plot_df, aes(x = pos, y = -value)) +
    geom_point(aes(colour = factor(codon_i), size = is_high),
               show.legend = c(color = TRUE, size = TRUE)) +
    scale_color_manual(name = "Codon position",
                       values = c("1" = "red", "2" = "orange", "3" = "green")) +
    scale_size_manual(name = "Key positions",
                     values = c(`FALSE` = 1.5, `TRUE` = 2.7)) +
    labs(x = "Nucleotide position", y = "Log-likelihood",
         title = sprintf("-log-likelihood across positions %d-%d (%s)",
                        min(index_rows), max(index_rows), gene_name)) +
    theme_minimal() + theme(legend.position = "right", legend.box = "vertical")
}

#' Calculate probability of correct base prediction

#' @param gene Gene variant name
#' @param df_gene Dataframe with predictions
#' @param all_variants DNAStringSet with reference sequences
#' @return String with probability result
probability_correct_base <- function(gene, df_gene, all_variants) {
  split_seq <- strsplit(as.character(all_variants[gene]), "")[[1]]
  correct_count <- 0
  
  for (i in 1:length(split_seq)) {
    max_column <- which.max(df_gene[i, 1:4])
    predicted_base <- c("A", "C", "G", "T")[max_column]
    if (predicted_base == split_seq[i]) {
      correct_count <- correct_count + 1  
    }
  }
  answer <- correct_count/length(split_seq)
  sprintf("The probability that EVO2 predicts the correct base across the gyrase gene in e-coli for %s is %.4f",
          gene, answer)
}

#' Plot metric differences between variants
#' 
#' @param ref_gene Reference gene name
#' @param target_gene Target gene name
#' @param ref_df Reference gene dataframe
#' @param target_df Target gene dataframe
#' @param metric Metric to compare ("entropy" or "log_likelihood")
#' @param index_rows Range of positions to plot
#' @param highlight Positions to highlight
#' @return ggplot object
plot_metric_diff <- function(ref_gene, target_gene, ref_df, target_df, 
                           metric = "entropy", index_rows = 200:300, highlight = 247:249, reverse_comp = FALSE) {
  if (reverse_comp){
    highlight <- (nrow(ref_df) - highlight + 1)
    index_rows <- sort(nrow(ref_df) - index_rows + 1)
  }
  diff_vec <- target_df[[metric]] - ref_df[[metric]]
  df <- data.frame(
    pos = index_rows,
    diff_vec = diff_vec[index_rows],
    codon_i = ((index_rows - 1) %% 3) + 1,
    is_high = index_rows %in% highlight
  )
  
  ggplot(df, aes(x = pos, y = diff_vec)) +
    geom_point(aes(colour = factor(codon_i), size = is_high),
               show.legend = c(colour = TRUE, size = TRUE)) +
    scale_color_manual("Codon position",
                      values = c("1" = "red", "2" = "orange", "3" = "green")) +
    scale_size_manual("Key positions",
                     values = c(`FALSE` = 1.5, `TRUE` = 2.7)) +
    labs(x = "Nucleotide position",
         y = sprintf("Δ %s", gsub("_", " ", metric, fixed = TRUE)),
         title = sprintf("%s difference: %s vs %s",
                        gsub("_", " ", metric, fixed = TRUE),
                        target_gene, ref_gene)) +
    theme_minimal() + theme(legend.position = "right", legend.box = "vertical")
}

#' Save plot as PDF
#' 
#' @param plot_obj ggplot object to save
#' @param filename Name for the PDF file
#' @param out_dir Output directory
#' @param width PDF width in inches
#' @param height PDF height in inches
save_as_pdf <- function(plot_obj, filename, out_dir = "My_figures", width = 12, height = 7) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  full_path <- file.path(out_dir, filename)
  
  ggsave(full_path,
         plot = plot_obj,
         device = cairo_pdf,   
         width = width,
         height = height,
         units = "in")
  message("Saved: ", full_path)
}


