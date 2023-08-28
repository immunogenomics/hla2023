#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

## This script will run a pseudobulk eQTL analysis on multiple datasets at once, conditioning on the lead, secondary, and tertiary variants.
## Input: PEER-normalized residuals for each dataset
## Output: HLA conditional (tertiary) eQTLs

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
    library(stringr)
})

message('read inputs')
cell_type = args[1] # Any of B_plasma, T, or Myeloid
# Read in residuals
resid = readRDS(paste0('/data/srlab1/jkang/hla2023/eqtl_pseudobulk/4_multidataset_eQTLs/', cell_type, '_residuals.rds'))
outprefix = paste0('/data/srlab1/jkang/hla2023/eqtl_pseudobulk/6_multidataset_conditional/', cell_type)
HLA_genes = c('HLA.A', 'HLA.B', 'HLA.C', 'HLA.DRB1', 'HLA.DQA1', 'HLA.DQB1', 'HLA.DPA1', 'HLA.DPB1')

# Reorder factor levels so that Randolph_NI is base level
resid$dataset = factor(resid$dataset, levels = c('Randolph_NI', 'AMP2RA', 'OneK1K', 'Smillie'))
print(paste('dimensions of samples X residuals matrix:', nrow(resid), 'X', ncol(resid)))

# Read in variants
geno_dosage_AMP2RA = readRDS('/data/srlab1/jkang/hla2023/data/sampleXdosage/AMP2RA_sampleXdosage_final.rds')
geno_dosage_OneK1K = readRDS('/data/srlab1/jkang/hla2023/data/sampleXdosage/OneK1K_sampleXdosage_final.rds')
geno_dosage_Smillie = readRDS('/data/srlab1/jkang/hla2023/data/sampleXdosage/Smillie2019_sampleXdosage_final.rds')
geno_dosage_Randolph_NI = readRDS('/data/srlab1/jkang/hla2023/data/sampleXdosage/Randolph_NI_sampleXdosage_final.rds')

# Combine all datasets into a single samples X dosage matrix
geno_dosage = rbind(geno_dosage_AMP2RA, geno_dosage_OneK1K, geno_dosage_Smillie, geno_dosage_Randolph_NI)
print(paste('dimensions of samples X dosages matrix:', nrow(geno_dosage), 'X', ncol(geno_dosage)))

# Add variants to the data frame
data = merge(resid, geno_dosage, by = 0)
data$Row.names = NULL # remove extra column added from call to merge

# Get lead variants for the cell type
lead_variants = read.csv(paste0('/data/srlab1/jkang/hla2023/eqtl_pseudobulk/4_multidataset_eQTLs/', 
                               cell_type, '_lead_variants.csv'))

secondary_lead_variants = read.csv(paste0(outprefix, '_secondary_lead_variants.csv'))

tertiary_lead_variants = read.csv(paste0(outprefix, '_tertiary_lead_variants.csv'))

# Run model
message('Start eQTL modeling')
full_results = NULL
mod = NULL
for (gene in HLA_genes) {
    message('Calculating eQTLs for gene: ', gene)
    
    # Get lead variants for the gene (to condition on)
    lead_variant = lead_variants[which(lead_variants$gene == str_replace(gene, '\\.', '-')), 'variant'] %>% unlist()
    secondary_lead_variant = secondary_lead_variants[which(secondary_lead_variants$gene == gene), 'variant'] %>% unlist()
    tertiary_lead_variant = tertiary_lead_variants[which(tertiary_lead_variants$gene == gene), 'variant'] %>% unlist()
    
    message('Conditioning on lead variant: ', lead_variant, secondary_lead_variant, tertiary_lead_variant)
    
    variants_to_test = colnames(data)[10:ncol(data)] # variants start in 10th row of data
    variants_to_test =  variants_to_test[variants_to_test != lead_variant &
                                         variants_to_test != secondary_lead_variant &
                                         variants_to_test != tertiary_lead_variant ] # remove leads
    
    results = rbindlist(lapply(variants_to_test, function(variant) { 
        mod = lm(as.formula(paste(gene, ' ~ dataset + ', lead_variant, '+', secondary_lead_variant, '+', tertiary_lead_variant, '+',
        paste(variant, collapse = "+"), sep = "")), data = data)
        
        # if variant is in perfect LD with a lead, then may not be possible to estimate beta
        if (variant %in% rownames(summary(mod)$coefficients)) {
            data.table(variant = variant,
                lead_variant = lead_variant,
                secondary_lead_variant = secondary_lead_variant,
                tertiary_lead_variant = tertiary_lead_variant,
                conditional_iter = 3,
                cell_type = cell_type,
                gene = gene,
                beta = summary(mod)$coefficients[variant, 'Estimate'], 
                stderr = summary(mod)$coefficients[variant, 'Std. Error'],
                t.val = summary(mod)$coefficients[variant, 't value'],
                p.val = summary(mod)$coefficients[variant, 'Pr(>|t|)'])
        }
    }))
    full_results = rbind(full_results, results)
}

message('Save results')
geno_df = readRDS('/data/srlab1/jkang/hla2023/data/sampleXdosage/four_cohorts_variant_info_final.rds')
geno_df$POS = as.numeric(geno_df$POS)
full_results = merge(full_results, geno_df, by.x = 'variant', by.y = 'ID', all.x = T, all.y = F)
lead_variants = full_results %>% group_by(gene) %>% slice(which.min(p.val))

write.csv(full_results, paste0(outprefix, '_quaternary_all_variants.csv'))
write.csv(lead_variants, paste0(outprefix, '_quaternary_lead_variants.csv'))