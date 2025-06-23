# inializing
library(Biostrings)
setwd("~/Desktop/Relman_lab/Intro_assignment")

source("example.r")

dna_set <- readDNAStringSet("gene_variants.fasta")
seqs <- as.character(dna_set)
seq_S1 <- seqs["83_S1"]
# done initializing


initialize_df <- function(EVO2_npy_file, sequence_name) { # should be .npy file and both inputs in ""
  prob_matrix <- get_logits(EVO2_npy_file) # normalizes the data into probabilities and outputs
  # a df with probabilities of nucleotides for each pos
  prob_matrix$entropy <- apply(prob_matrix[,1:4], 1, # calculates entropy and adds it as a column to the df
                               function(position_prob) -sum(position_prob * log2(position_prob)))
  # Need to get real FASTA sequences for log_likelihood
  dna_set <- readDNAStringSet("gene_variants.fasta") # makes a color coded object with (set of strings)
  # length + sequences + names out of FASTA file
  seqs <- as.character(dna_set)   # makes plain text vector and preserves name and sequence of genes
  df_genes <- data.frame(variant = # makes a data frame with the first column seq names 
                           # stringsAFactors = False tells R to keep strings (seqs) as strings 
                           names(dna_set),stringsAsFactors = FALSE)
  df_genes$sequence <- seqs # adds the sequence column to the genes df
  desired_seq <- seqs[sequence_name] # gets the desired sequence of the gene
  seq_vec <- strsplit(desired_seq, split = "")[[1]] # turns the sequence into a list of single nucleotides
  log_likelihood <- numeric(length(seq_vec)) # intializes a vector with as many 0's as length of the seq
  for (i in 1:length(seq_vec)) { # iterates from 1 to length of the sequence
    base <- seq_vec[i] # stores the correct base at position i
    base_index <- match(base, c("A", "C", "G", "T")) # matches the correct base to a column # in the df
    log_likelihood[i] <- log(prob_matrix[i, base_index]) # calculates the log likelihood at row i
  }
  prob_matrix$log_likelihood <- log_likelihood # adds the log likelihood column to prob matrix df
  
  return(prob_matrix)
}

# creates data frame with calculated entropy and log likelihood for S1 & other genes
df_S1 <- initialize_df("input_83_S1_logits.npy", "83_S1") # save data from S1 prob matrix
df_L <- initialize_df("input_83_L_logits.npy", "83_L") 
df_S2 <- initialize_df("input_83_S2_logits.npy", "83_S2") 
df_Stop <- initialize_df("input_83_Stop_logits.npy", "83_Stop") 

# Plot Entropy Function (positions 200-300 by default)
plot_entropy <- function(gene_name, gene_df, start = 200, end = 300, col = "blue") {
  positions <- start:end
  plot(positions, gene_df$entropy[positions],
       type = "l", col = col,
       xlab = "Nucleotide Position",
       ylab = "Entropy",
       main = sprintf("Entropy across positions %d–%d (%s)", start, end, gene_name))
}
plot_entropy("83_S1", df_S1, col = "green") # Plotting Entropy for "83_S1"

# plotting all 4 entropy plots for easy comparison
all_dfs <- list( # puts the df's and genes in a list w/ names = genes, values = df's 
  "83_S1"  = df_S1,
  "83_S2"  = df_S2,
  "83_L"   = df_L,
  "83_Stop"= df_Stop
)

cols <- c("blue", "darkgreen", "purple", "red") # creates a vector of colors

par(mfrow = c(2,2))  # put all four plots in a 2×2 grid
i <- 1 # to interate through color columns
for (name in names(all_dfs)) { # iterates through names of genes in all_dfs
  plot_entropy(name, all_dfs[[name]], col = cols[i]) # all_dfs[[name]] gives corresponding value
  i <- i + 1 # to iterate through colors
}
par(mfrow = c(1,1))  # reset plotting layout




# Plotting log-likelihood across positions 200-300 on 83_S1
positions <- 200:300 # using for all gene sequences
plot(positions, df_S1$log_likelihood[positions],
     type = "l", col = "blue",
     xlab = "Nucleotide Position",
     ylab = "Log-Likelihood",
     main = "Log-Likelihood across positions 200–300 (83_S1)")

# Calculating percentage of positions that most likely nucleotide (EVO2) was observed nucleotide (S1)

probability_correct_base <- function(gene, df_gene, dna_set, seq_S1){ # input is a gene string such as "83_S1" and corresponding df 
  correct_count <- 0
  split_seq <- strsplit(as.character(dna_set[gene]), "")[[1]]
  for (i in 1:length(split_seq)){
    max_column <- which.max(df_gene[i, 1:4])
    predicted_base <- c("A", "C", "G", "T")[max_column]
    if (predicted_base == seq_S1[i]){
      correct_count <- correct_count + 1  
    }
  }
  return(correct_count/length(split_seq))
}

probability_correct_base("83_S1", df_S1, dna_set, seq_S1)


# Initializing df for other gene variations
df_L <- initialize_df("input_83_L_logits.npy", "83_L") 
df_S2 <- initialize_df("input_83_S2_logits.npy", "83_S2") 
df_Stop <- initialize_df("input_83_Stop_logits.npy", "83_Stop") 


# Plotting Entropy across positions 200-300 on 83_L
plot(positions, df_L$entropy[positions],
     type = "l", col = "blue",
     xlab = "Nucleotide Position",
     ylab = "Entropy",
     main = "Entropy across positions 200–300 (83_L)")

# Plotting log-likelihood across positions 200-300 on 83_L
plot(positions, df_L$log_likelihood[positions],
     type = "l", col = "blue",
     xlab = "Nucleotide Position",
     ylab = "Log-Likelihood",
     main = "Log-Likelihood across positions 200–300 (83_L)")



# Plotting Entropy across positions 200-300 on 83_S2
plot(positions, df_S2$entropy[positions],
     type = "l", col = "blue",
     xlab = "Nucleotide Position",
     ylab = "Entropy",
     main = "Entropy across positions 200–300 (83_S2)")

# Plotting log-likelihood across positions 200-300 on 83_S2
plot(positions, df_S2$log_likelihood[positions],
     type = "l", col = "blue",
     xlab = "Nucleotide Position",
     ylab = "Log-Likelihood",
     main = "Log-Likelihood across positions 200–300 (83_S2)")



# Plotting Entropy across positions 200-300 on 83_Stop
plot(positions, df_Stop$entropy[positions],
     type = "l", col = "blue",
     xlab = "Nucleotide Position",
     ylab = "Entropy",
     main = "Entropy across positions 200–300 (83_Stop)")

# Plotting log-likelihood across positions 200-300 on 83_Stop
plot(positions, df_Stop$log_likelihood[positions],
     type = "l", col = "blue",
     xlab = "Nucleotide Position",
     ylab = "Log-Likelihood",
     main = "Log-Likelihood across positions 200–300 (83_Stop)")

# Now plotting per-basepair differences in log likelihood and entropy vs 83_S1
entropy_diff_L <- df_L$entropy - df_S1$entropy
loglike_diff_L <- df_L$log_likelihood - df_S1$log_likelihood

entropy_diff_S2 <- df_S2$entropy - df_S1$entropy
loglike_diff_S2 <- df_S2$log_likelihood - df_S1$log_likelihood

entropy_diff_Stop <- df_Stop$entropy - df_S1$entropy
loglike_diff_Stop <- df_Stop$log_likelihood - df_S1$log_likelihood

# Entropy diff plot 83_L
plot(positions, entropy_diff_L[positions],
     type = "l", col = "purple",
     xlab = "Nucleotide Position",
     ylab = "Δ Entropy",
     main = "Entropy Difference: 83_L vs 83_S1")

# Log-likelihood diff plot 83_L
plot(positions, loglike_diff_L[positions],
     type = "l", col = "darkgreen",
     xlab = "Nucleotide Position",
     ylab = "Δ Log-Likelihood",
     main = "Log-Likelihood Difference: 83_L vs 83_S1")


# Entropy diff plot 83_S2
plot(positions, entropy_diff_S2[positions],
     type = "l", col = "purple",
     xlab = "Nucleotide Position",
     ylab = "Δ Entropy",
     main = "Entropy Difference: 83_S2 vs 83_S1")

# Log-likelihood diff plot 83_S2
plot(positions, loglike_diff_S2[positions],
     type = "l", col = "darkgreen",
     xlab = "Nucleotide Position",
     ylab = "Δ Log-Likelihood",
     main = "Log-Likelihood Difference: 83_S2 vs 83_S1")


# Entropy diff plot 83_Stop
plot(positions, entropy_diff_Stop[positions],
     type = "l", col = "purple",
     xlab = "Nucleotide Position",
     ylab = "Δ Entropy",
     main = "Entropy Difference: 83_Stop vs 83_S1")

# Log-likelihood diff plot 83_Stop
plot(positions, loglike_diff_Stop[positions],
     type = "l", col = "darkgreen",
     xlab = "Nucleotide Position",
     ylab = "Δ Log-Likelihood",
     main = "Log-Likelihood Difference: 83_Stop vs 83_S1")

# Making 4x4 matrix of log-likelihood sums
total_ll <- c( #storing the sums in a vector
  "83_S1" = sum(df_S1$log_likelihood),
  "83_S2" = sum(df_S2$log_likelihood),
  "83_L"  = sum(df_L$log_likelihood),
  "83_Stop" = sum(df_Stop$log_likelihood)
)

diff_matrix <- outer(total_ll, total_ll, "-") # making a 4x4 matrix 
# comparing log likelihood
heatmap(diff_matrix, Rowv = NA, Colv = NA, scale = "none",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        main = "Log-Likelihood Differences Between Variants")












