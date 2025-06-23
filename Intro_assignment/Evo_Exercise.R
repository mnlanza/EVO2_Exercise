source("example.r")

setwd("~/Desktop/Relman lab/Intro assignment")

prob_matrix <- get_logits("input_83_S1_logits.npy")

#head(prob_matrix)# head() prints out first 6 rows

prob_matrix$entropy <- apply(prob_matrix, 1, 
              function(position_prob) -sum(position_prob * log2(position_prob)))
              # uses entropy formula to add a column of entropy for each position

'nrow(prob_matrix) # how many rows in matrix
cat(prob_matrix[200:300, 5], sep = "\n")# prints entropy for rows 200-300

max_entropy_val <- max(prob_matrix[200:300, 5]) # gets max entropy value
min_entropy_val <- min(prob_matrix[200:300, 5]) # gets min entropy value
max_entropy_row <- which.max(prob_matrix[200:300, 5]) # gets max entropy row (ordered from 1-101)
min_entropy_row <- which.min(prob_matrix[200:300, 5]) # gets min entropy row (ordered from 1-101)

cat("Max entropy value:", max_entropy_val,"at row", max_entropy_row + 199, "\n")
prob_matrix[max_entropy_row +199,]
cat("Min entropy value:", min_entropy_val, "at row", min_entropy_row + 199, "\n")
prob_matrix[min_entropy_row +199,]
'

      # Now calculating log likelihood with FASTA data

# install if needed:
if (!requireNamespace("Biostrings", quietly=TRUE))
  BiocManager::install("Biostrings")

library(Biostrings)

# Read the FASTA into a DNAStringSet
fasta_path <- "gene_variants.fasta"
dna_set <- readDNAStringSet(fasta_path)

# dna_set is a named vector of sequences
# names(dna_set) gives you the headers
# dna_set[[i]] is an XString object you can coerce to character()

names(dna_set)

#named character vector
seqs <- as.character(dna_set)
codons <- substr(seqs, 247, 249)
codons

df_genes <- data.frame( # does this data frame help?
  variant  = names(dna_set),
  stringsAsFactors = FALSE
)

df_genes$sequence <- seqs # appends the sequences to the data frame

seq_S1 <- seqs["83_S1"]
seq_vec <- strsplit(seq_S1, split = "")[[1]] # all of the nucleotides in 83_S1(FASTA) 
                                             # split into a list

log_likelihood <- numeric(length(seq_vec))

for (i in 1:length(seq_vec)) {
  base <- seq_vec[i]
  base_index <- match(base, c("A", "C", "G", "T"))
  log_likelihood[i] <- log(prob_matrix[i, base_index])
}

prob_matrix <- cbind(prob_matrix, log_likelihood)


correct_count <- 0
for (i in 1:length(seq_vec)){
  max_column <- which.max(prob_matrix[i, 1:4])
  predicted_base <- c("A", "C", "G", "T")[max_column]
  if (predicted_base == seq_vec[i]){
    correct_count <- correct_count + 1  
  }
}
probability_correct_base <- correct_count/length(seq_vec)

probability_correct_base
