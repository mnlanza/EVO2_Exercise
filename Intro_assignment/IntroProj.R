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

df_genes <- data.frame(
  variant  = names(dna_set),
  stringsAsFactors = FALSE
)

source("example.r")

setwd("~/Desktop/Relman lab/Intro assignment")


pp <- get_logits("input_83_S1_logits.npy")

head(pp)

pp$entropy <- apply(pp, 1, function(p) -sum(p * log2(p)))

nrow(pp)
pp[1:2,]
pp[500:502,]
pp[1000:1002,]
pp[1500:1502,]
pp[2000:2002,]
pp[2500:2502,]
pp[2623:2625,]

