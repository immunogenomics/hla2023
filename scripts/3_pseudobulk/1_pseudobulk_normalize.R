#!/usr/bin/Rscript

# Inputs:
# - cell type specific expression matrix (genes x cells) (.rds)
# - cell type specific metadata (cells x attributes) (.csv)
# - sample metadata (samples x attributes) (.csv)
# - prefix to store output

# Performs:
# log(CP10k+1) normalization within each cell
# Aggregates all cells per sample by mean log(CP10k+1) expression of each gene to get samplesXgenes matrix
# Rank-based inverse normal transforms each gene
# Output is ready for PEER

args = commandArgs(trailingOnly = TRUE)

message('loading libraries')
suppressPackageStartupMessages({
    library(tidyr)
    library(dplyr)
    library(plyr)
    library(Matrix)
    library(Matrix.utils)
    library(singlecellmethods)
})

exp = readRDS(args[1]) # path to cell type specific expression matrix (genes x cells)
cell_meta = readRDS(args[2]) # path to cell type specific metadata (cells x attributes)
sample_meta = read.csv(args[3], row.names = 1) # path to sample metadata (samples x attributes)
outdir = args[4] # prefix to store output

# Make pseudobulk profiles by single-cell log-normalization then mean aggregation across cells
exp_norm = normalizeData(exp, method = "log")
exp_norm = exp_norm[which(rowSums(exp_norm) != 0), ] # remove genes that are 0 across all cells
pseudobulk_scnorm_exp_sum = aggregate.Matrix(t(exp_norm), as.factor(cell_meta$Sample), fun = 'sum') # samples x genes
pseudobulk_scnorm_exp_mean = pseudobulk_scnorm_exp_sum / count(cell_meta$Sample)$freq

# Include only genes with non-zero expression in > half of the samples
pseudobulk_scnorm_exp_mean = pseudobulk_scnorm_exp_mean[, colSums(pseudobulk_scnorm_exp_mean > 0) > .5 * nrow(pseudobulk_scnorm_exp_mean)]
    
# Some samples might not have that cell type, so subset sample_meta
sample_meta_subset = sample_meta[rownames(pseudobulk_scnorm_exp_mean), ]
write.csv(sample_meta_subset, paste0(outdir, '_samples_meta.csv'), row.names = TRUE, quote = FALSE) # Write sample meta
if(all(rownames(pseudobulk_scnorm_exp_mean) != rownames(sample_meta_subset))) {
    message('Error: ordering of samples inconsistent.')}

# Write normalized result
saveRDS(pseudobulk_scnorm_exp_mean, paste0(outdir, '_samplesXgenes_norm.rds'))

# Rank-based inverse normal transformation
exp_int = apply(pseudobulk_scnorm_exp_mean, 2, function(x) {qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)) )})

# Write inverse normal transformed result
write.csv(exp_int, paste0(outdir, '_samplesXgenes_norm_invnt.csv'), row.names = TRUE, quote = FALSE)