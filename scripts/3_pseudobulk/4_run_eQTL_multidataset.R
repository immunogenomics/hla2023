#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

## Runs a pseudobulk eQTL analysis on multiple datasets at once using scHLApers expression estimates.
## Input: PEER-normalized residuals for each dataset
## Output: HLA eQTLs

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
})

cell_type = args[1] # Any of B_plasma, T, or Myeloid
prefix = '/data/srlab1/jkang/hla2023/'
outprefix = paste0(prefix, 'eqtl_pseudobulk/4_multidataset_eQTLs/', cell_type)

message('read pers inputs')
AMP2RA_resid_path = paste0(prefix, 'eqtl_pseudobulk/2_PEER/pers_AMP2RA_', cell_type, '_K7_PEER_residuals.csv')
OneK1K_resid_path = paste0(prefix,'eqtl_pseudobulk/2_PEER/pers_OneK1K_', cell_type, '_K20_PEER_residuals.csv')
Smillie_resid_path = paste0(prefix,'eqtl_pseudobulk/2_PEER/pers_Smillie_', cell_type, '_K2_PEER_residuals.csv')
Randolph_NI_resid_path = paste0(prefix,'eqtl_pseudobulk/2_PEER/pers_Randolph_NI_', cell_type, '_K7_PEER_residuals.csv')

AMP2RA_sXd_path = paste0(prefix, 'data/sampleXdosage/AMP2RA_sampleXdosage_final.rds')
OneK1K_sXd_path = paste0(prefix, 'data/sampleXdosage/OneK1K_sampleXdosage_final.rds')
Smillie_sXd_path = paste0(prefix, 'data/sampleXdosage/Smillie2019_sampleXdosage_final.rds')
Randolph_NI_sXd_path = paste0(prefix, 'data/sampleXdosage/Randolph_NI_sampleXdosage_final.rds')

# Read in residuals
AMP2RA_resid = read.csv(AMP2RA_resid_path, header = TRUE, sep = ',', row.names = 1) %>% as.matrix() %>% as.data.frame()
OneK1K_resid = read.csv(OneK1K_resid_path, header = TRUE, sep = ',', row.names = 1) %>% as.matrix() %>% as.data.frame()
Smillie_resid = read.csv(Smillie_resid_path, header = TRUE, sep = ',', row.names = 1) %>% as.matrix() %>% as.data.frame()
Randolph_NI_resid = read.csv(Randolph_NI_resid_path, header = TRUE, sep = ',', row.names = 1) %>% as.matrix() %>% as.data.frame()

# Subset residuals to HLA genes
HLA_genes = c('HLA.A', 'HLA.B', 'HLA.C', 'HLA.DRB1', 'HLA.DQA1', 'HLA.DQB1', 'HLA.DPA1', 'HLA.DPB1')
AMP2RA_resid = AMP2RA_resid[, HLA_genes]
OneK1K_resid = OneK1K_resid[, HLA_genes]
Smillie_resid = Smillie_resid[, HLA_genes]
Randolph_NI_resid = Randolph_NI_resid[, HLA_genes]

# Add a column to the residuals dataframe with dataset info
AMP2RA_resid$dataset = 'AMP2RA'
OneK1K_resid$dataset = 'OneK1K'
Smillie_resid$dataset = 'Smillie'
Randolph_NI_resid$dataset = 'Randolph_NI'

# Read in variants
geno_dosage_AMP2RA = readRDS(AMP2RA_sXd_path)[rownames(AMP2RA_resid), ] # sync sample names (not all cell types have all samples)
geno_dosage_OneK1K = readRDS(OneK1K_sXd_path)[rownames(OneK1K_resid), ]
geno_dosage_Smillie = readRDS(Smillie_sXd_path)[rownames(Smillie_resid), ]
geno_dosage_Randolph_NI = readRDS(Randolph_NI_sXd_path)[rownames(Randolph_NI_resid), ]

message('Combine datasets')
# Combine all datasets into a single samples X residuals matrix
resid = rbind(AMP2RA_resid, OneK1K_resid, Smillie_resid, Randolph_NI_resid)
print(paste('dimensions of samples X residuals matrix:', nrow(resid), 'X', ncol(resid)))

# Reorder factor levels so that Randolph_NI is base level
resid$dataset = factor(resid$dataset, levels = c('Randolph_NI', 'AMP2RA', 'OneK1K', 'Smillie'))

# Combine all datasets into a single samples X dosage matrix
geno_dosage = rbind(geno_dosage_AMP2RA, geno_dosage_OneK1K, geno_dosage_Smillie, geno_dosage_Randolph_NI)
#geno_dosage = geno_dosage[, which(colSums(geno_dosage) > 0)] ## Remove variants that do not have differences across individuals
print(paste('dimensions of samples X dosages matrix:', nrow(geno_dosage), 'X', ncol(geno_dosage)))

# Add variants to the data frame
data = merge(resid, geno_dosage, by = 0)
data$Row.names = NULL # remove extra column added from call to merge

# Run model
message('Start eQTL modeling')

full_results = NULL
for (gene in HLA_genes) {
    message('Calculating eQTLs for gene: ', gene)
    results = rbindlist(lapply(colnames(data)[10:ncol(data)], function(variant) { # variants start in 10th row of data
        mod = lm(as.formula(paste(gene, ' ~ dataset + ',
            paste(variant, collapse = "+"), sep = "")), data = data)
        
        data.table(variant = variant,
               cell_type = cell_type,
               gene = gene,
               beta = summary(mod)$coefficients[variant, 'Estimate'], 
               stderr = summary(mod)$coefficients[variant, 'Std. Error'],
               t.val = summary(mod)$coefficients[variant, 't value'],
               p.val = summary(mod)$coefficients[variant, 'Pr(>|t|)'])
    }))
    full_results = rbind(full_results, results)
}

message('Save results')
write.csv(full_results, paste0(outprefix, '_pseudobulk_eQTLs.csv'), quote = F)
saveRDS(resid, paste0(outprefix, '_residuals.rds'))
saveRDS(geno_dosage, paste0(outprefix, '_sampleXdosage.rds'))