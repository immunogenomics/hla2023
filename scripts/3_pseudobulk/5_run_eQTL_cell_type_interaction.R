#!/usr/bin/Rscript

## Script that will run a pseudobulk eQTL analysis on multiple datasets at once. This analysis will
## test for interaction term between genotype and cell type.
## Input: PEER-normalized residuals for each dataset
## Output: Significance of cell type interaction (across myeloid, B, and T cells)

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

message('read inputs')
outprefix = paste0('/data/srlab1/jkang/hla2023/eqtl_pseudobulk/5_multidataset_interaction/')

##########################################
#  Preprocess residuals for each cell type
##########################################
for (cell_type in c('B_plasma', 'Myeloid', 'T')) {
    message(cell_type)
    # Read in residuals
    resid = readRDS(paste0('/data/srlab1/jkang/hla2023/eqtl_pseudobulk/4_multidataset_eQTLs/', cell_type, '_residuals.rds'))
    
    # Reorder factor levels so that Randolph_NI is base level
    resid$dataset = factor(resid$dataset, levels = c('Randolph_NI', 'AMP2RA', 'OneK1K', 'Smillie'))
    
    print(paste('dimensions of samples X residuals matrix:', nrow(resid), 'X', ncol(resid)))
    resid$cell_type = cell_type
    saveRDS(resid, paste0(outprefix, cell_type, '_HLA_resid.rds'))
}

# Combine all cell types together into a single samples X residuals matrix
B_plasma_resid = readRDS(paste0(outprefix, 'B_plasma', '_HLA_resid.rds'))
B_plasma_resid$donor = rownames(B_plasma_resid)
Myeloid_resid = readRDS(paste0(outprefix, 'Myeloid', '_HLA_resid.rds'))
Myeloid_resid$donor = rownames(Myeloid_resid)
T_resid = readRDS(paste0(outprefix, 'T', '_HLA_resid.rds'))
T_resid$donor = rownames(T_resid)
all_resid = rbind(B_plasma_resid, Myeloid_resid, T_resid)

# Read in variants
geno_dosage_AMP2RA = readRDS('/data/srlab1/jkang/hla2023/data/sampleXdosage/AMP2RA_sampleXdosage_final.rds')
geno_dosage_OneK1K = readRDS('/data/srlab1/jkang/hla2023/data/sampleXdosage/OneK1K_sampleXdosage_final.rds')
geno_dosage_Smillie = readRDS('/data/srlab1/jkang/hla2023/data/sampleXdosage/Smillie2019_sampleXdosage_final.rds')
geno_dosage_Randolph_NI = readRDS('/data/srlab1/jkang/hla2023/data/sampleXdosage/Randolph_NI_sampleXdosage_final.rds')

# Combine all datasets into a single samples X dosage matrix
geno_dosage = rbind(geno_dosage_AMP2RA, geno_dosage_OneK1K, geno_dosage_Smillie, geno_dosage_Randolph_NI)
print(paste('dimensions of samples X dosages matrix:', nrow(geno_dosage), 'X', ncol(geno_dosage)))

# Add variants to the data frame
data = merge(all_resid, geno_dosage, by.x = 'donor', by.y = 0, all.x = TRUE) # left join to add dosage to resid
rownames(data) = data$Row.names
data$Row.names = NULL # remove extra column added from call to merge
data$cell_type = factor(data$cell_type, levels = c('T', 'B_plasma', 'Myeloid'))

# Run model
message('Start eQTL modeling')

# Get list of variants that were significant in the original analysis
T_lead = read.csv('/data/srlab1/jkang/hla2023/eqtl_pseudobulk/4_multidataset_eQTLs/T_lead_variants.csv') %>% mutate(cell_type = 'T')
B_plasma_lead = read.csv('/data/srlab1/jkang/hla2023/eqtl_pseudobulk/4_multidataset_eQTLs/B_plasma_lead_variants.csv')
Myeloid_lead = read.csv('/data/srlab1/jkang/hla2023/eqtl_pseudobulk/4_multidataset_eQTLs/Myeloid_lead_variants.csv')
all_lead = rbind(T_lead, B_plasma_lead, Myeloid_lead)

all_lead$lrt_pval = Inf
for (i in 1:nrow(all_lead)) {
    gene = str_replace(all_lead[i, 'gene'], '-', '.')
    variant = all_lead[i, 'variant']
    interaction_mod = lmer(as.formula(paste0(gene, ' ~ dataset + ', variant, " + ", 
                                                variant, "*cell_type + cell_type", "+ (1 | donor)")), data = data)
        
    null_mod = lmer(as.formula(paste0(gene, ' ~ dataset + ', variant, " + ", 
                                                "cell_type", "+ (1 | donor)")), data = data) # no interaction term
        
    all_lead$lrt_pval[i] = lrtest(interaction_mod, null_mod)$`Pr(>Chisq)`[2]
}
all_lead$lrt_qval = p.adjust(all_lead$lrt_pval, method = 'BH')

message('Save results')
write.csv(all_lead, paste0(outprefix, 'cell_type_interaction_lrt.csv'))