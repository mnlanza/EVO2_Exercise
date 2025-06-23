library(reticulate)

setwd("~/Desktop/Relman_lab/Intro_assignment") # set path

np <- import("numpy")

softmax <- function(x) { # takes a vector of logits (real #'s) and transforms into prob distr
  exp_x <- exp(x - max(x)) # subtracting each element from its max (for numerical stability)
  exp_x / sum(exp_x) # normalizing vector to sum to 1
}

get_logits <- function(fn) {
    base.np <- np$load(fn) # loading the function and putting into base.np

    # convert numpy file to Rx
    mat <- py_to_r(base.np)[1, , ] # python to r array (takes first slice)
    
    # turn into probability matrix using softmax (N samples x M tokens)
    all <- t(apply(mat, 1, softmax)) # applies 3rd arg to 1st arg (1 = row wise)
                                     # t is transpose
    # focus only on nucleotides 
    data.frame(A=all[,utf8ToInt("A") + 1], # utf8ToInt ASCII + of "A" and then adding 1
               C=all[,utf8ToInt("C") + 1],
               G=all[,utf8ToInt("G") + 1],
               T=all[,utf8ToInt("T") + 1])
}

# raw logits matrix (N samples x M tokens)
pp <- get_logits("input_83_S1_logits.npy")
pp
