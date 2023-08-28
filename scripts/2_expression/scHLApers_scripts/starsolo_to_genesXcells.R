#!/usr/bin/env Rscript
# Converts STARsolo output for 1 sample into a genesXcells matrix, subset to only include cells
# with barcodes matching the metadata "Cell" column.
# Inputs: 
#    - sample name
#    - results directory from STARsolo run
#    - string label added to prefix for run name (used for internal testing, can be removed)
#    - path to cell metadata

print('loading libraries')
suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(stats)
    library(tidyverse)
    library(Matrix)
    library(tidyr)
})

args = commandArgs(trailingOnly=TRUE)
sample = args[1] # sample name
STARsolo_results_dir = args[2] #e.g. '/data/srlab2/jkang/scHLA/personalized_final/Randolph2021_NewPanel/STARsolo_results/'
label=args[3] #e.g. "GeneFull_Ex50pAS"
cell_meta_path = args[4] # e.g. '/data/srlab1/jkang/hla/data/combined_AMP_Randolph_Smillie/cell_meta_Randolph.csv'
STARsolo_sample_dir = paste0(STARsolo_results_dir, sample, '_', label, '/')

print('read in list of good cells')
cell_meta = read.csv(cell_meta_path)

print('get cell and gene names for sample')
# get cell and gene names
result_cells = readLines(paste0(STARsolo_sample_dir, sample, '_', label, '_Solo.out/GeneFull_Ex50pAS/raw/barcodes.tsv'))
result_genes = readLines(paste0(STARsolo_sample_dir, sample, '_', label, '_Solo.out/GeneFull_Ex50pAS/raw/features.tsv'))

get_gene_name = function(gene) { return(str_split(gene, '\t')[[1]][2]) }
result_gene_names = unlist(lapply(result_genes, get_gene_name))

print('build starting matrices')
# build starting matrices
result_EM_mat = as(readMM(paste0(STARsolo_sample_dir, sample, '_', label, '_Solo.out/GeneFull_Ex50pAS/raw/UniqueAndMult-EM.mtx')), 'dgCMatrix') # EM 
dimnames(result_EM_mat) = list(result_gene_names, result_cells)
colnames(result_EM_mat) = paste0(sample, '_', colnames(result_EM_mat))

result_noMM_mat = as(readMM(paste0(STARsolo_sample_dir, sample, '_', label, '_Solo.out/GeneFull_Ex50pAS/raw/matrix.mtx')), 'dgCMatrix') # No multimapping
dimnames(result_noMM_mat) = list(result_gene_names, result_cells)
colnames(result_noMM_mat) = paste0(sample, '_', colnames(result_noMM_mat))

cells_sample = cell_meta$Cell[which(cell_meta$Sample == sample)] # Get list of cells for that sample

print('subset by QC cells')
# subset matrix by QC'd cells
result_EM_mat = result_EM_mat[, which(colnames(result_EM_mat) %in% cells_sample)]
result_noMM_mat = result_noMM_mat[, which(colnames(result_noMM_mat) %in% cells_sample)]

print('save!')
# save
saveRDS(result_EM_mat, paste0(STARsolo_sample_dir, 'exp_EM.rds'))
saveRDS(result_noMM_mat, paste0(STARsolo_sample_dir, 'exp_noMM.rds'))
