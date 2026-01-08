#!/usr/bin/env Rscript

# ==========================================================
# Cohesion calculation (aligns with your current algorithm)
# - persistence filter (pers.cutoff)
# - rel abundance
# - cor.mat.true
# - null model: for each focal taxon, permute others; take median expected cor
# - obs-exp correlations (or custom correlation matrix if use.custom.cors=TRUE)
# - connectedness: mean(pos), mean(neg)
# - cohesion: relab %*% connectedness
#
# Input: otu_all.txt (rows=samples, cols=taxa)
# Output:
#   out/connectedness_taxon.csv
#   out/cohesion_sample.csv
#
# Usage:
#   Rscript scripts/cohesion_calc.R --otu otu_all.txt --outdir out \
#     --pers_cutoff 0.10 --iter 200 --tax_shuffle TRUE \
#     --use_custom_cors FALSE
#
# ==========================================================

suppressPackageStartupMessages({
  library(optparse)
})

# -----------------------------
# parameters (CLI)
# -----------------------------
option_list <- list(
  make_option(c("--otu"), type="character", default="otu_all.txt",
              help="OTU table path (tab-delimited, header=TRUE, rownames=1)."),
  make_option(c("--outdir"), type="character", default="out",
              help="Output directory."),
  make_option(c("--pers_cutoff"), type="double", default=0.10,
              help="Persistence cutoff (fraction of samples with >0)."),
  make_option(c("--iter"), type="integer", default=200,
              help="Null iterations (>=200 recommended)."),
  make_option(c("--tax_shuffle"), type="character", default="TRUE",
              help="TRUE: shuffle each taxon across samples, keep focal fixed. FALSE: shuffle within-sample nonzero taxa (excluding focal)."),
  make_option(c("--use_custom_cors"), type="character", default="FALSE",
              help="TRUE to use custom correlation matrix (bypass null model subtraction)."),
  make_option(c("--custom_cor"), type="character", default=NULL,
              help="Custom correlation matrix csv path (header=TRUE, rownames=1)."),
  make_option(c("--seed"), type="integer", default=1,
              help="Random seed for reproducibility.")
)

opt <- parse_args(OptionParser(option_list=option_list))

pers.cutoff <- opt$pers_cutoff
iter <- opt$iter
tax.shuffle <- toupper(opt$tax_shuffle) == "TRUE"
use.custom.cors <- toupper(opt$use_custom_cors) == "TRUE"
set.seed(opt$seed)

if(!file.exists(opt$otu)) stop("OTU file not found: ", opt$otu)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# helper functions (same as yours)
# -----------------------------
zero <- function(vec){
  length(which(vec == 0))
}

neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals) == 0) n.mean <- 0
  n.mean
}

pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals) == 0) p.mean <- 0
  p.mean
}

# -----------------------------
# read data
# -----------------------------
b <- read.table(opt$otu, header = TRUE, row.names = 1, quote = "", sep = "\t", check.names = FALSE)
c <- as.matrix(b)
storage.mode(c) <- "numeric"

# remove all-zero samples/taxa
c <- c[rowSums(c) > 0, colSums(c) > 0, drop = FALSE]

# persistence cutoff in "number of samples present"
zero.cutoff <- ceiling(pers.cutoff * nrow(c))

# keep taxa present in at least zero.cutoff samples:
# apply(c,2,zero) counts zeros; require zeros < (n - zero.cutoff)
keep_taxa <- apply(c, 2, zero) < (nrow(c) - zero.cutoff)
d <- c[, keep_taxa, drop = FALSE]

# remove all-zero samples after filtering
d <- d[rowSums(d) > 0, , drop = FALSE]

message("Samples kept: ", nrow(d))
message("Taxa kept: ", ncol(d))

# -----------------------------
# relative abundance (FIXED: correct denominator)
# -----------------------------
rowsums.d <- rowSums(d)
rel.d <- d / rowsums.d

# -----------------------------
# custom correlation (optional)
# -----------------------------
if(use.custom.cors){
  if(is.null(opt$custom_cor) || !file.exists(opt$custom_cor)){
    stop("use_custom_cors=TRUE but --custom_cor is missing or not found.")
  }
  custom.cor.mat <- read.csv(opt$custom_cor, header = TRUE, row.names = 1, check.names = FALSE)
  custom.cor.mat <- as.matrix(custom.cor.mat)
  
  # subset to retained taxa and align order
  taxa <- colnames(rel.d)
  if(!all(taxa %in% colnames(custom.cor.mat))){
    stop("Custom correlation matrix does not contain all retained taxa.")
  }
  custom.cor.mat.sub <- custom.cor.mat[taxa, taxa, drop = FALSE]
  
  # use custom correlations directly as obs-exp matrix (same as your ifelse logic)
  obs.exp.cors.mat <- custom.cor.mat.sub
  diag(obs.exp.cors.mat) <- 0
  
} else {
  
  # -----------------------------
  # true correlation
  # -----------------------------
  cor.mat.true <- cor(rel.d, use = "pairwise.complete.obs")
  diag(cor.mat.true) <- 0
  
  # -----------------------------
  # null model (MEDIAN expected cors per focal taxon)
  # -----------------------------
  med.tax.cors <- NULL
  
  if(tax.shuffle){
    # tax.shuffle=TRUE: shuffle each taxon across samples, keep focal fixed
    for(which.taxon in 1:ncol(rel.d)){
      perm.cor.vec.mat <- NULL
      
      for(i in 1:iter){
        # create permuted matrix
        perm.rel.d <- matrix(0, nrow(rel.d), ncol(rel.d))
        rownames(perm.rel.d) <- rownames(rel.d)
        colnames(perm.rel.d) <- colnames(rel.d)
        
        # permute every taxon across samples
        for(j in 1:ncol(rel.d)){
          perm.rel.d[, j] <- sample(rel.d[, j])
        }
        # keep focal taxon unshuffled
        perm.rel.d[, which.taxon] <- rel.d[, which.taxon]
        
        cor.mat.null <- cor(perm.rel.d, use = "pairwise.complete.obs")
        diag(cor.mat.null) <- 0
        
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
      }
      
      # median expected correlation between focal and all taxa
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median, na.rm = TRUE))
      
      if(which.taxon %% 20 == 0) message("Processed focal taxon: ", which.taxon)
    }
    
  } else {
    # tax.shuffle=FALSE: shuffle within each sample among nonzero taxa excluding focal
    for(which.taxon in 1:ncol(rel.d)){
      perm.cor.vec.mat <- NULL
      
      for(i in 1:iter){
        perm.rel.d <- rel.d
        
        for(j in 1:nrow(rel.d)){
          which.replace <- which(rel.d[j, ] > 0)
          which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]
          if(length(which.replace.nonfocal) > 1){
            perm.rel.d[j, which.replace.nonfocal] <- sample(rel.d[j, which.replace.nonfocal])
          }
        }
        
        cor.mat.null <- cor(perm.rel.d, use = "pairwise.complete.obs")
        diag(cor.mat.null) <- 0
        
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
      }
      
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median, na.rm = TRUE))
      
      if(which.taxon %% 20 == 0) message("Processed focal taxon: ", which.taxon)
    }
  }
  
  colnames(med.tax.cors) <- colnames(rel.d)
  rownames(med.tax.cors) <- colnames(rel.d)
  
  # -----------------------------
  # obs - expected
  # -----------------------------
  obs.exp.cors.mat <- cor.mat.true - med.tax.cors
  diag(obs.exp.cors.mat) <- 0
}

# -----------------------------
# connectedness (taxon-level)
# -----------------------------
connectedness.pos <- apply(obs.exp.cors.mat, 2, pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat, 2, neg.mean)

connectedness_df <- data.frame(
  taxon = names(connectedness.pos),
  Positive_Connectedness = as.numeric(connectedness.pos),
  Negative_Connectedness = as.numeric(connectedness.neg),
  stringsAsFactors = FALSE
)

# -----------------------------
# cohesion (sample-level)
# -----------------------------
cohesion.pos <- as.numeric(rel.d %*% connectedness.pos)
cohesion.neg <- as.numeric(rel.d %*% connectedness.neg)

cohesion_df <- data.frame(
  sample = rownames(rel.d),
  Positive_Cohesion = cohesion.pos,
  Negative_Cohesion = cohesion.neg,
  stringsAsFactors = FALSE
)

# -----------------------------
# write outputs (FIXED: split files)
# -----------------------------
write.csv(connectedness_df, file = file.path(opt$outdir, "connectedness_taxon.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(cohesion_df, file = file.path(opt$outdir, "cohesion_sample.csv"),
          row.names = FALSE, quote = FALSE)

message("Done. Wrote:\n- ", file.path(opt$outdir, "connectedness_taxon.csv"),
        "\n- ", file.path(opt$outdir, "cohesion_sample.csv"))
