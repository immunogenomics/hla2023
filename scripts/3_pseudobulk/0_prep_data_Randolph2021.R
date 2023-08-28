#!/usr/bin/Rscript

# Perform cell QC for Randolph datasets
# And write cell type specific matrices

# Inputs:
# - Gene expression matrix (genes x cells) (.rds)
# - Cell metadata (cells x attributes) (.csv)
# - sample metadata (samples x attributes) (.csv)
# - prefix to store output

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

prefix = '/data/srlab1/jkang/hla2023/data/'
prefix_out = '/data/srlab1/jkang/hla2023/eqtl_pseudobulk/'

###################
## Randolph2021 ###
###################
message('Randolph2021')
# This cohort is special due to the repeated samples per individual. 
# To process this sample, let's split by condition (flu/NI) into two "datasets".

message('Read in data')
cell_meta = read.csv(paste0(prefix, 'meta/cell_meta_Randolph_completeHLA.csv'), row.names = 1)
genesXcells_pers = readRDS(paste0(prefix, 'genesXcells_GeneFull_Ex50pAS/exp_Randolph2021_pers_EM_GeneFull_Exon50pAS.rds'))

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

idx_keep = which(cell_meta$percent.mito < 0.2 & cell_meta$nGene >= 500)
cell_meta = cell_meta[idx_keep, ]

nremoved = (ncol(genesXcells_pers) - length(idx_keep))
print(paste('Cells before filtering:', ncol(genesXcells_pers), 'cells'))
print(paste('Removed', nremoved, 'cells'))
print(paste('Cells remaining:', nrow(cell_meta)))

# Save filtered cell meta file with QC metrics
write.csv(cell_meta, paste0(prefix, 'meta/cell_meta_Randolph_completeHLA_cellQC.csv'), quote = F)


# Read in data
sample_meta = read.csv(paste0(prefix, 'meta/sample_meta_Randolph_completeHLA.csv'), row.names = 1)
cell_meta = read.csv(paste0(prefix, 'meta/cell_meta_Randolph_completeHLA_cellQC.csv'), row.names = 1)
genesXcells_pers = readRDS(paste0(prefix, 'genesXcells_GeneFull_Ex50pAS/exp_Randolph2021_pers_EM_GeneFull_Exon50pAS.rds'))
genesXcells_noPers = readRDS(paste0(prefix, 'genesXcells_GeneFull_Ex50pAS/exp_Randolph2021_noPers_noMM_GeneFull_Exon50pAS.rds'))

# Rename 'IMGT' to 'HLA'
rownames(genesXcells_pers) = str_replace(rownames(genesXcells_pers), 'IMGT_', 'HLA-')

# Fix duplicated row (gene) names
rownames(genesXcells_pers)[which(duplicated(rownames(genesXcells_pers)))]
rownames(genesXcells_pers) = make.unique(rownames(genesXcells_pers)) # fix duplicate rownames
rownames(genesXcells_noPers)[which(duplicated(rownames(genesXcells_noPers)))]
rownames(genesXcells_noPers) = make.unique(rownames(genesXcells_noPers)) # fix duplicate rownames

## Subset cell_meta to the samples included in sample metadata
dim(cell_meta)
cell_meta = cell_meta %>% filter(Sample %in% sample_meta$Sample)
dim(cell_meta)

# Divide metadata into flu and NI
sample_meta_NI = sample_meta[which(endsWith(sample_meta$Sample, 'NI')), ]
sample_meta_flu = sample_meta[which(endsWith(sample_meta$Sample, 'flu')), ]
cell_meta_NI = cell_meta[which(endsWith(cell_meta$Sample, 'NI')), ]
cell_meta_flu = cell_meta[which(endsWith(cell_meta$Sample, 'flu')), ]
dim(cell_meta_NI)
dim(cell_meta_flu)

## Process NI

# Subset and match ordering of cells/samples
genesXcells_pers_NI = genesXcells_pers[, cell_meta_NI$Cell] # Subset and match the ordering of cells in genesXcells to cell_meta
genesXcells_noPers_NI = genesXcells_noPers[, cell_meta_NI$Cell] # Subset and match the ordering of cells in genesXcells to cell_meta
dim(genesXcells_pers_NI)
dim(genesXcells_noPers_NI)

# Write cell-type-specific matrices
table(cell_meta_NI$cell_type_major)
for (c in unique(cell_meta_NI$cell_type_major)) {
    message(c)
    idx_c = which(cell_meta_NI$cell_type_major == c)
    
    cell_meta_NI_c = cell_meta_NI[idx_c, ]
    genesXcells_pers_NI_c = genesXcells_pers_NI[, idx_c]
    genesXcells_noPers_NI_c = genesXcells_noPers_NI[, idx_c]
    
    # Remove individuals with fewer than 5 cells of this cell type
    cells_per_sample = cell_meta_NI_c %>% group_by(Sample) %>% dplyr::summarise(n = n())
    samples_less5cells = cells_per_sample$Sample[which(cells_per_sample$n < 5)]
    message('Sample has fewer than 5 cells: ', samples_less5cells)
    idx_keep = which(! cell_meta_NI_c$Sample %in% samples_less5cells)
    
    cell_meta_NI_c = cell_meta_NI_c[idx_keep, ]
    message('Num cells: ', nrow(cell_meta_NI_c))
    message('Num samples: ', length(unique(cell_meta_NI_c$Sample)))
    saveRDS(cell_meta_NI_c[idx_keep, ], paste0(prefix_out, '0_data/cell_meta_Randolph_NI_', c, '.rds'))
    saveRDS(genesXcells_pers_NI_c[, idx_keep], paste0(prefix_out, '0_data/exp_pers_Randolph_NI_', c, '.rds'))
    saveRDS(genesXcells_noPers_NI_c[, idx_keep], paste0(prefix_out, '0_data/exp_noPers_Randolph_NI_', c, '.rds'))
}

## Process flu

# Subset and match ordering of cells/samples
genesXcells_pers_flu = genesXcells_pers[, cell_meta_flu$Cell] # Subset and match the ordering of cells in genesXcells to cell_meta
genesXcells_noPers_flu = genesXcells_noPers[, cell_meta_flu$Cell] # Subset and match the ordering of cells in genesXcells to cell_meta
dim(genesXcells_pers_flu)
dim(genesXcells_noPers_flu)

# Write cell-type-specific matrices
table(cell_meta_flu$cell_type_major)
for (c in unique(cell_meta_flu$cell_type_major)) {
    message(c)
    idx_c = which(cell_meta_flu$cell_type_major == c)
    
    cell_meta_flu_c = cell_meta_flu[idx_c, ]
    genesXcells_pers_flu_c = genesXcells_pers_flu[, idx_c]
    genesXcells_noPers_flu_c = genesXcells_noPers_flu[, idx_c]
    
    # Remove individuals with fewer than 5 cells of this cell type
    cells_per_sample = cell_meta_flu_c %>% group_by(Sample) %>% dplyr::summarise(n = n())
    samples_less5cells = cells_per_sample$Sample[which(cells_per_sample$n < 5)]
    message('Sample has fewer than 5 cells: ', samples_less5cells)
    idx_keep = which(! cell_meta_flu_c$Sample %in% samples_less5cells)
    
    cell_meta_flu_c = cell_meta_flu_c[idx_keep, ]
    message('Num cells: ', nrow(cell_meta_flu_c))
    message('Num samples: ', length(unique(cell_meta_flu_c$Sample)))
    saveRDS(cell_meta_flu_c[idx_keep, ], paste0(prefix_out, '0_data/cell_meta_Randolph_flu_', c, '.rds'))
    saveRDS(genesXcells_pers_flu_c[, idx_keep], paste0(prefix_out, '0_data/exp_pers_Randolph_flu_', c, '.rds'))
    saveRDS(genesXcells_noPers_flu_c[, idx_keep], paste0(prefix_out, '0_data/exp_noPers_Randolph_flu_', c, '.rds'))
}

## Save NI and flu meta
write.csv(sample_meta_NI, paste0(prefix_out, '0_data/sample_meta_Randolph_NI.csv'), row.names = T, quote = F)
write.csv(sample_meta_flu, paste0(prefix_out, '0_data/sample_meta_Randolph_flu.csv'), row.names = T, quote = F)
write.csv(cell_meta_NI, paste0(prefix_out, '0_data/cell_meta_Randolph_NI.csv'), row.names = T, quote = F)
write.csv(cell_meta_flu, paste0(prefix_out, '0_data/cell_meta_Randolph_flu.csv'), row.names = T, quote = F)

message('All done!')