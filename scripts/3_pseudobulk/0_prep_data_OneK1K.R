#!/usr/bin/Rscript

# Perform cell QC for OneK1K dataset (too large for interactive session)
# And write cell type specific matrices
# Note: For all other cohorts, see 1_prep_data_submit_jobs.ipynb for cell QC

# Inputs:
# - Gene expression matrix (genes x cells) (.rds)
# - Cell metadata (cells x attributes) (.csv)
# - sample metadata (samples x attributes) (.csv)
# - prefix to store output

args = commandArgs(trailingOnly = TRUE)

message('loading libraries')
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggpubr)
    library(tidyr)
    library(dplyr)
    library(plyr)
    library(Matrix)
    library(vcfR)
    library(lme4)
    library(Matrix.utils)
    library(singlecellmethods)
    library(MASS)
    library(glmnet)
    library(patchwork)
    library(stringr)
    library(purrr)
    library(fitdistrplus)
    require(gdata)
    library(readxl)
})

fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width, repr.plot.res = 200)
}

mito.genes.coding = c('MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8','MT-ATP6','MT-CO3',
                      'MT-ND3','MT-ND4L','MT-ND4','MT-ND5','MT-ND6','MT-CYB')


message('Read in data')
cell_meta = read.csv('/data/srlab1/jkang/hla2023/data/meta/cell_meta_OneK1K_completeHLA.csv', row.names = 1)
genesXcells_pers = readRDS('/data/srlab1/jkang/hla2023/data/genesXcells_GeneFull_Ex50pAS/exp_OneK1K_pers_EM_GeneFull_Exon50pAS.rds')
out_gXc_dir = '/data/srlab1/jkang/hla2023/eqtl_pseudobulk/0_data/'

# Subset the genesXcells matrix to remove doublets
genesXcells_pers = genesXcells_pers[, cell_meta$Cell]

# Fix duplicate row (gene) names
rownames(genesXcells_pers)[which(duplicated(rownames(genesXcells_pers)))]
rownames(genesXcells_pers) = make.unique(rownames(genesXcells_pers)) 

message('Calculate nUMI and percent mitochondrial')
cell_meta$nUMI = colSums(genesXcells_pers)
cell_meta$percent.mito = Matrix::colSums(genesXcells_pers[mito.genes.coding, ]) / Matrix::colSums(genesXcells_pers)

message('Calculate nGene')
logical_count = (genesXcells_pers != 0)
cell_meta$nGene = Matrix::colSums(logical_count)
rm(logical_count)

# % MT, nUMI, nGene
print(paste('Max %MIT: ', max(cell_meta$percent.mit)))

print(paste('Min nUMI: ', round(min(cell_meta$nUMI))))
print(paste('Max nUMI: ', round(max(cell_meta$nUMI))))

print(paste('Min nGene: ', min(cell_meta$nGene)))
print(paste('Max nGene: ', max(cell_meta$nGene)))

# Perform filtering
idx_keep = which(cell_meta$percent.mito < 0.2 & cell_meta$nGene >= 500)
cell_meta = cell_meta[idx_keep, ]

nremoved = (ncol(genesXcells_pers) - length(idx_keep))
print(paste('Cells before filtering:', ncol(genesXcells_pers), 'cells'))
print(paste('Removed', nremoved, 'cells'))
print(paste('Cells remaining:', nrow(cell_meta)))

# Plot distributions for % MT, nUMI, nGene AFTER FILTERING
print('AFTER FILTERING')
print(paste('Max %MIT: ', max(cell_meta$percent.mit)))

print(paste('Min nUMI: ', round(min(cell_meta$nUMI))))
print(paste('Max nUMI: ', round(max(cell_meta$nUMI))))

print(paste('Min nGene: ', min(cell_meta$nGene)))
print(paste('Max nGene: ', max(cell_meta$nGene)))

# Save filtered cell meta file with QC metrics
write.csv(cell_meta, '/data/srlab1/jkang/hla2023/data/meta/cell_meta_OneK1K_completeHLA_cellQC.csv', quote = F)


message('Making cell type specific matrices')
# Read in data
sample_meta = read.csv('/data/srlab1/jkang/hla2023/data/meta/sample_meta_OneK1K_completeHLA.csv', row.names = 1)
cell_meta = read.csv('/data/srlab1/jkang/hla2023/data/meta/cell_meta_OneK1K_completeHLA_cellQC.csv', row.names = 1)
genesXcells_pers = readRDS('/data/srlab1/jkang/hla2023/data/genesXcells_GeneFull_Ex50pAS/exp_OneK1K_pers_EM_GeneFull_Exon50pAS.rds')
genesXcells_noPers = readRDS('/data/srlab1/jkang/hla2023/data/genesXcells_GeneFull_Ex50pAS/exp_OneK1K_noPers_noMM_GeneFull_Exon50pAS.rds')

## Subset cell_meta to the samples included in sample metadata
dim(cell_meta)
cell_meta = cell_meta %>% filter(Sample %in% sample_meta$Sample)
dim(cell_meta)

# Rename 'IMGT' to 'HLA'
rownames(genesXcells_pers) = str_replace(rownames(genesXcells_pers), 'IMGT_', 'HLA-')

# Fix duplicated row (gene) names
rownames(genesXcells_pers)[which(duplicated(rownames(genesXcells_pers)))]
rownames(genesXcells_pers) = make.unique(rownames(genesXcells_pers)) # fix duplicate rownames
rownames(genesXcells_noPers)[which(duplicated(rownames(genesXcells_noPers)))]
rownames(genesXcells_noPers) = make.unique(rownames(genesXcells_noPers)) # fix duplicate rownames

# Subset and match ordering of cells/samples
genesXcells_pers = genesXcells_pers[, cell_meta$Cell] # Subset and match the ordering of cells in genesXcells to cell_meta
genesXcells_noPers = genesXcells_noPers[, cell_meta$Cell] # Subset and match the ordering of cells in genesXcells to cell_meta
dim(genesXcells_pers)
dim(genesXcells_noPers)

# Write cell-type-specific matrices
table(cell_meta$cell_type_major)
for (c in unique(cell_meta$cell_type_major)) {
    message(c)
    idx_c = which(cell_meta$cell_type_major == c)
    
    cell_meta_c = cell_meta[idx_c, ]
    genesXcells_pers_c = genesXcells_pers[, idx_c]
    genesXcells_noPers_c = genesXcells_noPers[, idx_c]
    
    # Remove individuals with fewer than 5 cells of this cell type
    cells_per_sample = cell_meta_c %>% group_by(Sample) %>% dplyr::summarise(n = n())
    samples_less5cells = cells_per_sample$Sample[which(cells_per_sample$n < 5)]
    message('Sample has fewer than 5 cells: ', paste(samples_less5cells))
    idx_keep = which(! cell_meta_c$Sample %in% samples_less5cells)
    
    cell_meta_c = cell_meta_c[idx_keep, ]
    message('Num cells: ', nrow(cell_meta_c))
    message('Num samples: ', length(unique(cell_meta_c$Sample)))
    saveRDS(cell_meta_c, paste0(out_gXc_dir, 'cell_meta_OneK1K_', c, '.rds'))
    saveRDS(genesXcells_pers_c[, idx_keep], paste0(out_gXc_dir, 'exp_pers_OneK1K_', c, '.rds'))
    saveRDS(genesXcells_noPers_c[, idx_keep], paste0(out_gXc_dir, 'exp_noPers_OneK1K_', c, '.rds'))
}

message('Complete!')