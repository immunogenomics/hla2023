#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

## Script that will run a pseudobulk eQTL analysis on OneK1K samples. 
## This analysis will test for interaction term between genotype and age (Gxage) and sex (Gxsex)
## Input: 
#   - Cell type
#   - eGene
#   - SNP (eQTL variant to test for interaction)

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
})

# Arguments - test the lead SNPs for age interaction effects
cell_type = args[1]
gene = args[2] # HLA.DQA1
snp = args[3]

# Read in the sample-level expression and metadata
prefix = '/data/srlab1/jkang/hla2023/eqtl_pseudobulk/'
message('read inputs')
outprefix = paste0(prefix, '5_multidataset_interaction/', cell_type, '_', gene, '_', snp)

########################### OneK1K ################################
K = 20
# Read in OneK1K
exp_file = paste0(prefix, '1_normalization/pers_OneK1K_', cell_type, '_samplesXgenes_norm_invnt.csv')
meta_file = paste0(prefix, '1_normalization/pers_OneK1K_', cell_type, '_samples_meta.csv')
expr = read.csv(exp_file, header = TRUE, sep = ',', row.names = 1)
meta = read.csv(meta_file, header = TRUE, sep = ',', row.names = 1)

# Read in variants
sXd_path = paste0('/data/srlab1/jkang/hla2023/data/sampleXdosage/OneK1K_sampleXdosage_final.rds')
geno_dosage = readRDS(sXd_path)[rownames(expr), snp]

meta$snp = geno_dosage
meta$gene = expr[, gene]

# Regress out the PEER factors from the expression
peer_factors = read.csv(paste0(prefix, '2_PEER/pers_OneK1K_', cell_type, '_K', K, '_PEER_factors.csv'))
peer_factors = peer_factors[,(ncol(peer_factors) - K + 1):ncol(peer_factors)]
colnames(peer_factors) = paste0('PEER_', 1:20)
meta_OneK1K = cbind(meta, peer_factors)

mod = lm(as.formula(paste0('gene', ' ~ ',  paste0('PEER_', 1:K, collapse = '+'))), data = meta_OneK1K)
meta_OneK1K$exp_resid = resid(mod)

##################################################################

covs = c('Sample', 'Female', 'Age', 'gPC1', 'gPC2', 'gPC3', 'gPC4', 'gPC5', 'snp', 'exp_resid')
meta = meta_OneK1K[, covs]

message('Start eQTL modeling')

############### Test age and sex interaction
base_mod = lm(as.formula(paste0('exp_resid', ' ~ ',  
                                paste0(covs[2:9], collapse = '+'))), data = meta)
age_int_mod = lm(as.formula(paste0('exp_resid', ' ~ ',  'Age*snp + ',
                                paste0(covs[2:9], collapse = '+'))), data = meta)

sex_int_mod = lm(as.formula(paste0('exp_resid', ' ~ ',  'Female*snp + ',
                                paste0(covs[2:9], collapse = '+'))), data = meta)

res = list()

res[['base_mod']] = base_mod
res[['age_int_mod']] = age_int_mod
res[['age_lrt_P']] = lrtest(age_int_mod, base_mod)$`Pr(>Chisq)`[2]
res[['sex_int_mod']] = sex_int_mod
res[['sex_lrt_P']] = lrtest(sex_int_mod, base_mod)$`Pr(>Chisq)`[2]
res[['data']] = meta

message('Save results')
saveRDS(res, paste0(outprefix, '_age_sex_int_OneK1K.rds'))