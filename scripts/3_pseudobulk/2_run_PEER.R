#!/usr/bin/Rscript
# Runs PEER (Stegle et al. 2010)
# Use peer conda env

# Inputs: 
# - expression (samples x genes)
# - covariates (samples x metadata)
# - prefix for output files
# Performs peer normalization using K hidden factors

args = commandArgs(trailingOnly = TRUE)

message('loading libraries')
suppressPackageStartupMessages({
    library(peer)
})

message('read inputs')
exp_file = args[1]    # path to expression (samples x genes)
meta_file = args[2]   # path to covariates
out_prefix = args[3]  # output file prefix
K = args[4]           # number of hidden determinants to calculate

expr = read.csv(exp_file, header = TRUE, sep = ',', row.names = 1)
meta = read.csv(meta_file, header = TRUE, sep = ',', row.names = 1)
covs = cbind(meta[, 'Female'], meta[, 'Age'], 
             meta[, 'gPC1'], meta[, 'gPC2'], meta[, 'gPC3'], meta[, 'gPC4'], meta[, 'gPC5'])

colnames(covs) = c('female', 'age', 'gPC1', 'gPC2', 'gPC3', 'gPC4', 'gPC5')

# Add 10x chemistry batch variable if present
if ('Chemistry' %in% colnames(meta)) {
    covs = cbind(covs, model.matrix(~0+meta$Chemistry)[, 1])
    colnames(covs)[length(colnames(covs))] = 'chemistry'
}

message('set up PEER model')
model = PEER()
PEER_setPhenoMean(model, as.matrix(expr))
PEER_setAdd_mean(model, TRUE)
PEER_setNk(model, K)
PEER_getNk(model)
PEER_setCovariates(model, as.matrix(covs)) # Add covariates
PEER_setNmax_iterations(model, 10000)

message('Run PEER model')
PEER_update(model)

message('Save results')
residuals = PEER_getResiduals(model)
rownames(residuals) = rownames(expr)
colnames(residuals) = colnames(expr)
factors = PEER_getX(model)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)

# Save                        
write.csv(residuals, paste0(out_prefix, '_PEER_residuals.csv'), row.names = TRUE, quote = FALSE)
write.csv(factors, paste0(out_prefix, '_PEER_factors.csv'), row.names = TRUE, quote = FALSE)
write.csv(weights, paste0(out_prefix, '_PEER_weights.csv'), row.names = TRUE, quote = FALSE)
write.csv(precision, paste0(out_prefix, '_PEER_precision.csv'), row.names = TRUE, quote = FALSE)

# Plot convergence
pdf(paste0(out_prefix, "_PEER_plotModel.pdf"), width=8, height=8)
PEER_plotModel(model)
dev.off()

# Plot factor importance
Alpha = PEER_getAlpha(model)
write.csv(Alpha, paste0(out_prefix, '_PEER_alpha.csv'), row.names = TRUE, quote = FALSE)
pdf(paste0(out_prefix, "_PEER_alpha.pdf"), width=8, height=8)
plot(1.0 / Alpha, xlab = "Factors", ylab = "Factor relevance", main="")
dev.off()