#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

## Script that will run a pseudobulk eQTL analysis on OneK1K. This analysis will
## test for interaction term between eQTL and percent of cytotoxic / naive T cells.

message('loading libraries')
suppressPackageStartupMessages({
    library(data.table)
    library(tidyr)
    library(dplyr)
    library(plyr)
    library(Matrix)
    library(lme4)
    library(Matrix.utils)
    library(MASS)
    library(glmnet)
    require(gdata)
    library(lmtest)
    library(stringr)
})
source('/data/srlab1/jkang/hla2023/scripts/utils.R')

cell_type = 'T'
gene = args[1]
snp = args[2]

message('read inputs')
prefix = paste0('/data/srlab1/jkang/hla2023/eqtl_pseudobulk/')
outprefix = paste0(prefix, '5_multidataset_interaction/', cell_type, '_', gene, '_', snp)

# Read in residuals
resid = readRDS(paste0(prefix, '5_multidataset_interaction/', cell_type, '_HLA_resid.rds'))
resid$donor = rownames(resid)

# Read in variants
geno_dosage = readRDS(paste0(prefix, '/4_multidataset_eQTLs/', cell_type, '_sampleXdosage.rds'))

# Add variants to the data frame
data = merge(resid, geno_dosage, by.x = 'donor', by.y = 0, all.x = TRUE) # left join to add dosage to resid
rownames(data) = data$donor
colnames(data)[1] = 'Sample'

# Subset to OneK1K
data = data %>% filter(dataset == 'OneK1K')

# Read in percent cytotoxic cells
percent_clusters = readRDS(paste0(prefix, '/5_multidataset_interaction/T_cell_prop_cell_clusters.rds'))

data = left_join(data, percent_clusters, by = 'Sample', all.x = TRUE) # left join to add dosage to resid
data = data[, c('Sample', gene, snp, 'prop_cytotoxic', 'prop_naive')]


message('Start eQTL modeling')
############### Test interactions

res = list()

# Cytotoxic
res[['Cytotoxic_base_mod']] = lm(as.formula(paste0(gene, ' ~ prop_cytotoxic + ', snp)), data = data)
res[['Cytotoxic_int_mod']] = lm(as.formula(paste0(gene, ' ~ prop_cytotoxic + ', snp, 
                                                  ' + prop_cytotoxic*', snp)), data = data)
res[['Cytotoxic_lrt_P']] = lrtest(res[['Cytotoxic_base_mod']], res[['Cytotoxic_int_mod']])$`Pr(>Chisq)`[2]

# Naive
res[['Naive_base_mod']] = lm(as.formula(paste0(gene, ' ~ prop_naive + ', snp)), data = data)
res[['Naive_int_mod']] = lm(as.formula(paste0(gene, ' ~  prop_naive + ', snp, 
                                                  ' + prop_naive*', snp)), data = data)
res[['Naive_lrt_P']] = lrtest(res[['Naive_base_mod']], res[['Naive_int_mod']])$`Pr(>Chisq)`[2]

res[['data']] = data

message('Save results')
saveRDS(res, paste0(outprefix, '_T_cell_cluster_int.rds'))