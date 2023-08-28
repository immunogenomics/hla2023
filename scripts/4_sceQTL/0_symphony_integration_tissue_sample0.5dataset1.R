#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

# This script will create a Symphony object from multiple datasets for a major cell type by first integrating the two tissue datasets using Harmony, then mapping the two PBMC datasets onto the tissue-defined embedding using Symphony.
# Harmonizes over sample and dataset for the reference datasets, then maps the query each dataset at a time, by correcting for sample.

# Load libraries
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

# Inputs and constants
cell_type = args[1] # major cell type
prefix = '/data/srlab1/jkang/hla2023/'
outprefix = paste0(prefix, 'symphony/Tissue_', cell_type, '_sample0.5dataset1') # file names label

HLA_genes = c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 
              'CIITA', 'NLRC5', 'B2M', 'HLA-DMB', 'HLA-DMA', 'CD74', 'HLA-DRA', 'ACTB')

exp_combined = NULL
meta_combined = NULL

message('Read in data')
# Read in counts matrix
exp_AMP2RA = readRDS(paste0(prefix, 'eqtl_pseudobulk/0_data/exp_pers_AMP2RA_', cell_type, '.rds'))
exp_Smillie = readRDS(paste0(prefix, 'eqtl_pseudobulk/0_data/exp_pers_Smillie_', cell_type, '.rds'))

# Read in cell metadata
meta_AMP2RA = readRDS(paste0(prefix, 'eqtl_pseudobulk/0_data/cell_meta_AMP2RA_', cell_type, '.rds'))
meta_Smillie = readRDS(paste0(prefix, 'eqtl_pseudobulk/0_data/cell_meta_Smillie_', cell_type, '.rds'))

### Read in sample metadata to add Batch column
# Smillie - 10x chemistry
sample_meta_Smillie = read.csv(paste0(prefix, 'data/meta/sample_meta_Smillie_completeHLA.csv'))
sample_meta_Smillie$Batch = paste0(sample_meta_Smillie$Dataset, '_', sample_meta_Smillie$Chemistry) 
meta_Smillie = left_join(meta_Smillie, sample_meta_Smillie[, c('Sample', 'Batch')]) # Accidentally replaced rownames with numbers
rownames(meta_Smillie) = meta_Smillie$Cell
    
# AMP2RA - sample as batch
meta_AMP2RA$Batch = meta_AMP2RA$Sample
    
message('Combine datasets')
# Combine all datasets into a single expression matrix
exp_combined = cbind(exp_AMP2RA, exp_Smillie)
common_cols = intersect(colnames(meta_AMP2RA), colnames(meta_Smillie))
meta_combined = rbind(meta_AMP2RA[, common_cols], meta_Smillie[, common_cols])

dim(exp_combined)
dim(meta_combined)

message('Normalize data')
# Standard normalization and dimensionality reduction pipeline
exp_norm = singlecellmethods::normalizeData(exp_combined, 10000, 'log')

message('Save expression matrics')
saveRDS(exp_combined, paste0(prefix, 'symphony/Tissue_', cell_type, '_exp.rds'))
saveRDS(exp_norm, paste0(prefix, 'symphony/Tissue_', cell_type, '_exp_norm.rds'))

message('Variable gene selection')
ngenes = 1500
var_genes = singlecellmethods::vargenes_vst(exp_combined, groups = meta_combined$dataset, topn = ngenes)

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
        data_mat = t(Z_pca_ref),             ## PCA embedding matrix of cells
        meta_data = meta_combined,           ## dataframe with cell labels
        theta = c(0.5, 1),                   ## cluster diversity enforcement
        vars_use = c('Sample', 'dataset'),   ## variable to integrate out
        nclust = 50,                         ## number of clusters in Harmony model **
        sigma = 0.2,                         ## soft cluster fuzziness **
        max.iter.harmony = 50,
        return_object = TRUE,                ## return the full Harmony model object
        do_pca = FALSE                       ## don't recompute PCs
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

####################################
# Map PBMC datasets into tissue axes
####################################

tissue_ref = readRDS(paste0(outprefix, '_reference.rds'))

# Map OneK1K
query_exp = readRDS(paste0(prefix, 'eqtl_pseudobulk/0_data/exp_pers_OneK1K_', cell_type, '.rds'))
query_meta = readRDS(paste0(prefix, 'eqtl_pseudobulk/0_data/cell_meta_OneK1K_', cell_type, '.rds'))
queryOneK1K = mapQuery(exp_query = query_exp, metadata_query = query_meta, ref_obj = tissue_ref, vars = 'Sample')

# Map Randolph2021
query_exp_NI = readRDS(paste0(prefix, 'eqtl_pseudobulk/0_data/exp_pers_Randolph_NI_', cell_type, '.rds'))
query_meta_NI = readRDS(paste0(prefix, 'eqtl_pseudobulk/0_data/cell_meta_Randolph_NI_', cell_type, '.rds'))
query_exp_flu = readRDS(paste0(prefix, 'eqtl_pseudobulk/0_data/exp_pers_Randolph_flu_', cell_type, '.rds'))
query_meta_flu = readRDS(paste0(prefix, 'eqtl_pseudobulk/0_data/cell_meta_Randolph_flu_', cell_type, '.rds'))
query_exp = cbind(query_exp_NI, query_exp_flu)
query_meta = rbind(query_meta_NI, query_meta_flu)

queryRandolph = mapQuery(exp_query = query_exp, metadata_query = query_meta, ref_obj = tissue_ref, vars = 'Sample')
queryRandolph$meta_data = queryRandolph$meta_data %>% mutate(dataset = 
                                        ifelse(dataset %in% c('Randolph'), 'Randolph2021', dataset))

# Make ref/query object for plotting function
common_cols = intersect(colnames(tissue_ref$meta_data), colnames(queryOneK1K$meta_data))
ref_query = list()
ref_query$meta_data = rbind(tissue_ref$meta_data[, common_cols], 
                                queryOneK1K$meta_data[, common_cols], queryRandolph$meta_data[, common_cols])
ref_query$umap$embedding = rbind(tissue_ref$umap$embedding, queryOneK1K$umap, queryRandolph$umap)

ref_query$meta_data$dataset_cell_type_fine = paste0(ref_query$meta_data$dataset, '_', ref_query$meta_data$cell_type_fine)
unique(ref_query$meta_data$dataset_cell_type_fine)
# Group Randolph datasets together for faceting
ref_query$meta_data$dataset = factor(ref_query$meta_data$dataset, 
                                     levels = c('AMP2RA', 'Smillie2019', 'Randolph2021', 'OneK1K')) # reorder factor levels

exp_OneK1K = readRDS(paste0(prefix, 'symphony/OneK1K_', cell_type, '_batch2_reference.rds'))$exp_norm_HLA
exp_Randolph = readRDS(paste0(prefix, 'symphony/Randolph_', cell_type, '_batch2_reference.rds'))$exp_norm_HLA
ref_query$exp_norm_HLA = cbind(tissue_ref$exp_norm_HLA, exp_OneK1K, exp_Randolph)

# Save the ref/query object
ref_query$Z_corr = cbind(tissue_ref$Z_corr, queryOneK1K$Z, queryRandolph$Z)
saveRDS(ref_query, paste0(prefix, 'symphony/RefQuery_', cell_type,'_mapBloodOntoTissue.rds'))