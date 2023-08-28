#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

# Runs pseudobulk eQTL analysis on a single dataset

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

message('read inputs')
resid_path = args[1]   # Path to PEER residuals (.csv) 
sXd_path = args[2]     # Path to sampleXdosage file (.rds) 
outprefix = args[3]    # Prefix for output files

# Read in residuals
resid = read.csv(resid_path, header = TRUE, sep = ',', row.names = 1) %>% as.matrix() %>% as.data.frame()

# Subset residuals to HLA genes
HLA_genes = c('HLA.A', 'HLA.B', 'HLA.C', 'HLA.DRB1', 'HLA.DQA1', 'HLA.DQB1', 'HLA.DPA1', 'HLA.DPB1')
resid = resid[, HLA_genes]

# Read in variants
geno_dosage = readRDS(sXd_path)
geno_dosage = geno_dosage[rownames(resid), ] # sync sample names (not all cell types have all samples)

# Remove variants that do not have differences across individuals
geno_dosage = geno_dosage[, which(colSums(geno_dosage) > 0)]

# Add variants to the data frame
data = merge(resid, geno_dosage, by = 0)
data$Row.names = NULL # remove extra column added from call to merge

# Run eQTL model
message('Start eQTL modeling')
full_results = NULL
for (gene in HLA_genes) {
    message('Calculating eQTLs for gene: ', gene)
    results = rbindlist(lapply(colnames(data)[9:ncol(data)], function(variant) { 
        mod = lm(as.formula(paste(gene, ' ~ ', paste(variant, collapse = "+"), sep = "")), data = data)
        
        data.table(variant = variant, 
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