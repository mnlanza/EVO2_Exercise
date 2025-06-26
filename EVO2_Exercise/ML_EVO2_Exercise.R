# initializing
library(Biostrings)
setwd("~/Desktop/Relman_lab/EVO2_Exercise")

source("handle_EVO2_output.r")

dna_set <- readDNAStringSet("raw_data/gene_variants.fasta")
seqs <- as.character(dna_set)
seq_S1 <- seqs["83_S1"]
# done initializing


initialize_df <- function(EVO2_npy_file, sequence_name) { # should be .npy file and both inputs in ""
  prob_matrix <- get_logits(EVO2_npy_file) # normalizes the data into probabilities and outputs
  # a df with probabilities of nucleotides for each pos
  prob_matrix$entropy <- apply(prob_matrix[,1:4], 1, # calculates entropy and adds it as a column to the df
                               function(position_prob) -sum(position_prob * log2(position_prob)))
  # Need to get real FASTA sequences for log_likelihood
  dna_set <- readDNAStringSet("raw_data/gene_variants.fasta") # makes a color coded object with (set of strings)
  # length + sequences + names out of FASTA file
  seqs <- as.character(dna_set)   # makes plain text vector and preserves name and sequence of genes
  desired_seq <- seqs[sequence_name] # gets the desired sequence of the gene
  seq_vec <- strsplit(desired_seq, split = "")[[1]] # turns the sequence into a list of single nucleotides
  log_likelihood <- numeric(length(seq_vec)) # initializes a vector with as many 0's as length of the seq
  for (i in 1:length(seq_vec)) { # iterates from 1 to length of the sequence
    base <- seq_vec[i] # stores the correct base at position i
    base_index <- match(base, c("A", "C", "G", "T")) # matches the correct base to a column # in the df
    log_likelihood[i] <- log2(prob_matrix[i, base_index]) # calculates the log likelihood at row i
  }
  prob_matrix$log_likelihood <- log_likelihood # adds the log likelihood column to prob matrix df
  
  return(prob_matrix)
}

# creates data frame with calculated entropy and log likelihood for S1 & other genes
df_S1 <- initialize_df("raw_data/input_83_S1_logits.npy", "83_S1") # save data from S1 prob matrix
df_L <- initialize_df("raw_data/input_83_L_logits.npy", "83_L") 
df_S2 <- initialize_df("raw_data/input_83_S2_logits.npy", "83_S2") 
df_Stop <- initialize_df("raw_data/input_83_Stop_logits.npy", "83_Stop") 

library(ggplot2)

plot_entropy <- function(gene_name, gene_df, start = 200, end = 300, col = "blue") {
  EVO2_rows <- start:end  # storing start and end in positions
  ggplot(mapping = aes(x = start:end, y = gene_df$entropy[EVO2_rows])) +
    geom_point(color = col, size = 1.8) +
    labs(
      x = "Nucleotide Position",
      y = "Entropy",
      title = sprintf("Entropy across positions %d–%d (%s)", start, end, gene_name)
    ) +
    theme_minimal()
}
plot_entropy("83_S1", df_S1, col = "green") # Plotting Entropy for "83_S1"


library(patchwork) # for comparing entropy graphs

entropy_p1 <- plot_entropy("83_S1", df_S1, col = "blue") 
entropy_p2 <- plot_entropy("83_S2", df_S2, col = "darkgreen")
entropy_p3 <- plot_entropy("83_L",  df_L,  col = "purple")
entropy_p4 <- plot_entropy("83_Stop", df_Stop, col = "red")

(entropy_p1 | entropy_p2) / (entropy_p3 | entropy_p4)  # 2x2 layout using patchwork


# Plotting log-likelihood across positions 200-300 on 83_S1
library(ggplot2)
plot_log_likelihood <- function(gene_name,gene_df, start = 200, end = 300, col = "green"){
  EVO2_rows <- start:end  # account for pos 1 assoc w/ nuc 2 from EVO2 data
  
  ggplot(mapping = aes(x = start:end, y = gene_df$log_likelihood[EVO2_rows])) +
    geom_point(color = col, size = 1.8) +
    labs(
      x = "Nucleotide Position",
      y = "Log-likelihood",
      title = sprintf("Log-likelihood across positions %d–%d (%s)", start, end, gene_name)
    ) +
    theme_minimal()
}

# To compare log-likelihood graphs
log_ll_p1 <- plot_log_likelihood("83_S1", df_S1, col = "blue") 
log_ll_p2 <- plot_log_likelihood("83_S2", df_S2, col = "green") 
log_ll_p3 <- plot_log_likelihood("83_Stop", df_Stop, col = "red") 
log_ll_p4 <- plot_log_likelihood("83_L", df_L, col = "purple") 

(log_ll_p1 | log_ll_p2) / (log_ll_p3 | log_ll_p4)  # 2x2 layout using patchwork
  
# Calculating percentage of positions that most likely nucleotide (EVO2) was observed nucleotide (S1)

probability_correct_base <- function(gene, df_gene, dna_set){ # input is a gene string such as "83_S1" and corresponding df 
  correct_count <- 0
  split_seq <- strsplit(as.character(dna_set[gene]), "")[[1]]
  for (i in 1:length(split_seq)){
    max_column <- which.max(df_gene[i, 1:4])
    predicted_base <- c("A", "C", "G", "T")[max_column]
    if (predicted_base == split_seq[i]){
      correct_count <- correct_count + 1  
    }
  }
  return(correct_count/length(split_seq))
}

probability_correct_base("83_S1", df_S1, dna_set)
probability_correct_base("83_S2", df_S2, dna_set)
probability_correct_base("83_Stop", df_S1, dna_set)
probability_correct_base("83_L", df_S1, dna_set)


# Now plotting per-basepair differences in log likelihood and entropy vs 83_S1
entropy_diff_L <- df_L$entropy - df_S1$entropy
loglike_diff_L <- df_L$log_likelihood - df_S1$log_likelihood

entropy_diff_S2 <- df_S2$entropy - df_S1$entropy
loglike_diff_S2 <- df_S2$log_likelihood - df_S1$log_likelihood

entropy_diff_Stop <- df_Stop$entropy - df_S1$entropy
loglike_diff_Stop <- df_Stop$log_likelihood - df_S1$log_likelihood

positions <- (200:300) 
# Entropy diff plot 83_L (position 200-300)
plot(positions, entropy_diff_L[positions],
     type = "p", col = "purple",
     xlab = "Nucleotide Position",
     ylab = "Δ Entropy",
     main = "Entropy Difference: 83_L vs 83_S1",
     pch = 16)

# Log-likelihood diff plot 83_L (position 200-300)
plot(positions, loglike_diff_L[positions],
     type = "p", col = "darkgreen",
     xlab = "Nucleotide Position",
     ylab = "Δ Log-Likelihood",
     main = "Log-Likelihood Difference: 83_L vs 83_S1",
     pch = 16)


# Entropy diff plot 83_S2 (position 200-300)
plot(positions, entropy_diff_S2[positions],
     type = "p", col = "blue",
     xlab = "Nucleotide Position",
     ylab = "Δ Entropy",
     main = "Entropy Difference: 83_S2 vs 83_S1",
     pch = 16)

# Log-likelihood diff plot 83_S2 (position 200-300)
plot(positions, loglike_diff_S2[positions],
     type = "p", col = "darkgreen",
     xlab = "Nucleotide Position",
     ylab = "Δ Log-Likelihood",
     main = "Log-Likelihood Difference: 83_S2 vs 83_S1",
     pch = 16)


# Entropy diff plot 83_Stop (position 200-300)
plot(positions, entropy_diff_Stop[positions],
     type = "p", col = "purple",
     xlab = "Nucleotide Position",
     ylab = "Δ Entropy",
     main = "Entropy Difference: 83_Stop vs 83_S1",
     pch = 16)

# Log-likelihood diff plot 83_Stop (position 200-300)
plot(positions, loglike_diff_Stop[positions],
     type = "p", col = "darkgreen",
     xlab = "Nucleotide Position",
     ylab = "Δ Log-Likelihood",
     main = "Log-Likelihood Difference: 83_Stop vs 83_S1",
     pch = 16)


# Making 4x4 matrix of log-likelihood sums
total_ll <- c( #storing the sums in a vector
  "83_S1" = sum(df_S1$log_likelihood),
  "83_S2" = sum(df_S2$log_likelihood),
  "83_L"  = sum(df_L$log_likelihood),
  "83_Stop" = sum(df_Stop$log_likelihood)
)

diff_matrix <- outer(total_ll, total_ll, "-") # making a 4x4 matrix 


library(ggplot2)

df_plot <- as.data.frame(as.table(diff_matrix))   # rows: Var1, Var2, Freq

ggplot(df_plot, aes(Var1, Var2, fill = Freq)) + # Creates heatmap of log likelihood diff
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", Freq))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "Δ log-L") +
  labs(title = "Log-Likelihood Differences Between Variants",
       x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



