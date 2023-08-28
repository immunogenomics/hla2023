#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

# This script will create a Symphony object from the Randolph2021 cells.

suppressPackageStartupMessages({
    library(tidyr)
    library(dplyr)
    library(plyr)
    library(Matrix)
    library(Matrix.utils)
    library(singlecellmethods)
    library(irlba)
    library(symphony)
    library(ggplot2)
    library(RColorBrewer)
})

# Inputs
cell_type = args[1]
batch_theta = args[2]
outprefix = paste0('/data/srlab1/jkang/hla2023/symphony/Randolph_', cell_type, '_batch', batch_theta) # file names label

HLA_genes = c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 
              'CIITA', 'NLRC5', 'B2M', 'HLA-DMB', 'HLA-DMA', 'CD74', 'HLA-DRA', 'ACTB')
exp_combined = NULL
meta_combined = NULL

message('Read in data')
# Read in counts matrix
exp_flu = readRDS(paste0('/data/srlab1/jkang/hla2023/eqtl_pseudobulk/0_data/exp_pers_Randolph_flu_', cell_type, '.rds'))
exp_NI = readRDS(paste0('/data/srlab1/jkang/hla2023/eqtl_pseudobulk/0_data/exp_pers_Randolph_NI_', cell_type, '.rds'))
exp_combined = cbind(exp_flu, exp_NI)

# Read in cell metadata
meta_flu = readRDS(paste0('/data/srlab1/jkang/hla2023/eqtl_pseudobulk/0_data/cell_meta_Randolph_flu_', cell_type, '.rds'))
meta_NI = readRDS(paste0('/data/srlab1/jkang/hla2023/eqtl_pseudobulk/0_data/cell_meta_Randolph_NI_', cell_type, '.rds'))
sample_meta = read.csv('/data/srlab1/jkang/hla2023/data/meta/sample_meta_Randolph_completeHLA.csv')
meta_Randolph = rbind(meta_flu, meta_NI)
meta_combined = left_join(meta_Randolph, sample_meta[,-1])

dim(exp_combined)
dim(meta_combined)

message('Normalize data')
# Standard normalization and dimensionality reduction pipeline
exp_norm = singlecellmethods::normalizeData(exp_combined, 10000, 'log')

message('Variable gene selection')
var_genes = singlecellmethods::vargenes_vst(exp_combined, topn = 2000) # top 2000 vargenes

# Check if VST exclude cell cycle and MT and ribosome genes
mt_ri <- grep("^MT-|^RPL|^RPS|MALAT1|MIR-", row.names(exp_combined), value = TRUE)
cycle_prolif <- c(Seurat::cc.genes$s.genes, Seurat::cc.genes$g2m.genes)
genes_exclude <- union(mt_ri, cycle_prolif)
var_genes <- var_genes[-which(var_genes %in% genes_exclude)]
length(var_genes)

exp_vargenes = exp_norm[var_genes, ]

message('Start PCA')
# Z_score and PCA
vargenes_means_sds = tibble(symbol = var_genes, mean = Matrix::rowMeans(exp_vargenes))
vargenes_means_sds$stddev = singlecellmethods::rowSDs(exp_vargenes, vargenes_means_sds$mean)
head(vargenes_means_sds)
ref_exp_scaled = singlecellmethods::scaleDataWithStats(exp_vargenes, vargenes_means_sds$mean, vargenes_means_sds$stddev, 1)

set.seed(0)
s = irlba(ref_exp_scaled, nv = 10, maxit = 10000)
Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
loadings = s$u

message('Start Harmony')
# Harmony integration
set.seed(0)
ref_harmObj = harmony::HarmonyMatrix(
        data_mat = t(Z_pca_ref),           ## PCA embedding matrix of cells
        meta_data = meta_combined,         ## dataframe with cell labels
        theta = as.numeric(batch_theta),   ## cluster diversity enforcement
        vars_use = c('Sample'),            ## variable to integrate out
        nclust = 50,                       ## number of clusters in Harmony model
        sigma = 0.2,
        max.iter.harmony = 50,
        return_object = TRUE,              ## return the full Harmony model object
        do_pca = FALSE                     ## don't recompute PCs
)

message('Start Symphony')
# Compress a Harmony object into a Symphony reference
reference = symphony::buildReferenceFromHarmonyObj(
                           ref_harmObj,            # output object from HarmonyMatrix()
                           meta_combined,          # reference cell metadata
                           vargenes_means_sds,     # gene names, means, and std devs for scaling
                           loadings,               # genes x PCs matrix
                           verbose = TRUE,         # verbose output
                           do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
                           save_uwot_path = paste0(outprefix, '_uwot'),
                           umap_min_dist = 0.2)

message('Save HLA expression')
# Save HLA expression
reference$exp_HLA = exp_combined[HLA_genes, ]
reference$exp_norm_HLA = exp_norm[HLA_genes, ]

message('Make PCA UMAP')
# Make pre-Harmony UMAP
umap = uwot::umap(t(reference$Z_orig), n_neighbors = 30, learning_rate = 0.5, init = "laplacian", 
            metric = 'cosine', fast_sgd = FALSE, n_sgd_threads = 1, # for reproducibility
            min_dist = 0.2, n_threads = 4, ret_model = TRUE)
reference$pca_umap$embedding = umap$embedding
colnames(reference$pca_umap$embedding) = c('UMAP1', 'UMAP2')

message('Save final reference object')
# Save final reference
saveRDS(reference, paste0(outprefix, '_reference.rds'))

message('All done!')