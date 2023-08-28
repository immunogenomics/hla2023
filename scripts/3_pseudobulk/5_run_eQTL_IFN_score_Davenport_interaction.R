#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

## Script that will run score each sample by IFN signature, then perform a pseudobulk eQTL analysis to test for interaction term between genotype and IFN score.
## Input: 
#   - Cell type
#   - eGene
#   - SNP (eQTL variant to test for interaction)
#   - dataset ('OneK1K')
## Output: Significance of IFN-score interaction

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

cell_type = args[1]
gene = args[2] # HLA.DQA1
snp = args[3]
mydataset = args[4]

prefix = '/data/srlab1/jkang/hla2023/eqtl_pseudobulk/'
outprefix = paste0(prefix, '5_multidataset_interaction/', cell_type, '_', gene, '_', snp, '_', mydataset)

########### Score the samples
# Calculate signature for each sample
dataset_file = ifelse(mydataset == 'Randolph', 'Randolph_NI', mydataset)
exp_norm = readRDS(paste0(prefix, '1_normalization/pers_', dataset_file, 
                                 '_', cell_type, '_samplesXgenes_norm.rds'))

# Scale the sum of scaled expression (within each dataset)
score_IFN = scale(rowSums(scale(as.matrix(exp_norm[, which(colnames(exp_norm) %in% IFN_sig_Davenport)])))) %>% as.data.frame()
colnames(score_IFN) = 'IFN_score'
score_IFN$Sample = rownames(score_IFN)

########### Run the eQTL interaction
# Read in residuals
resid = readRDS(paste0(prefix, '5_multidataset_interaction/', cell_type, '_HLA_resid.rds')) %>%
    filter(dataset == dataset_file)
resid$donor = rownames(resid)

# Read in variants
geno_dosage = readRDS(paste0(prefix, '/4_multidataset_eQTLs/', cell_type, '_sampleXdosage.rds'))

# Add variants to the data frame
data = merge(resid, geno_dosage, by.x = 'donor', by.y = 0, all.x = TRUE) # left join to add dosage to resid
rownames(data) = data$donor
colnames(data)[1] = 'Sample'
data = merge(data, score_IFN, by = 'Sample', all.x = TRUE) # left join to add dosage to resid
data = data[, c('Sample', 'dataset', gene, snp, 'IFN_score')]


message('Start eQTL modeling')
############### Test interactions

res = list()
res[['IFN_score_base_mod']] = lm(as.formula(paste0(gene, ' ~ IFN_score + ', snp)), data = data)
res[['IFN_score_int_mod']] = lm(as.formula(paste0(gene, ' ~ IFN_score + ', snp, 
                                                  ' + IFN_score*', snp)), data = data)
res[['IFN_score_lrt_P']] = lrtest(res[['IFN_score_base_mod']], res[['IFN_score_int_mod']])$`Pr(>Chisq)`[2]
res[['data']] = data

message('Save results')
saveRDS(res, paste0(outprefix, '_IFN_score_Davenport_int.rds'))