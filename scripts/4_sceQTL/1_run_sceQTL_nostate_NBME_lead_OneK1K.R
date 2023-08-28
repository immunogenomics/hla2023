#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

# Given a cell type, raw counts, a single variant, and gene, this script will test whether the regulatory variant is an eQTL using the NBME model for the OneK1K dataset.
# Does not incorporate cell state (hPCs), as this tests for the genotype main effect only.

message('loading libraries')
suppressPackageStartupMessages({
    library(tidyr)
    library(dplyr)
    library(plyr)
    library(Matrix)
    library(lme4)
    library(Matrix.utils)
    library(singlecellmethods)
    library(MASS)
    library(glmnet)
    library(stringr)
    library(lmtest)
    library(lme4)
})

message('Read in single-cell data')
cell_type = args[1] #'T'
ref_path = paste0('/data/srlab1/jkang/hla2023/symphony/OneK1K_', cell_type, '_batch2_reference.rds') # path to ref object
gene = args[2] # 'HLA-DQA1'
snp = args[3] #'rs3104371'
ref = readRDS(ref_path)

# Load raw UMI counts matrix and metadata and gene expression PCs from Symphony object
meta = ref$meta_data
pca_res = ref$Z_orig        # pre-Harmony PCs
harmonypca_res = ref$Z_corr # post-Harmony PCs

message('Read in genotype dosages')
# Load genotype dosages
geno = readRDS('/data/srlab1/jkang/hla2023/data/sampleXdosage/OneK1K_sampleXdosage_final.rds')

message('Make dataframe')
# Make data frame of variables for model
E = as.numeric(ref$exp_HLA[gene, ]) %>% round() # round to get integers
G = t(geno[, snp])[, as.character(meta$Sample)]
IND = as.factor(meta$Geno_ID)
B = paste0('OneK1K_', meta$Batch)
AGE = scale(meta$Age)
SEX = meta$Female
nUMI = scale(log(meta$nUMI))
MT = scale(meta$percent.mito)

dataset = meta$dataset
expPC = pca_res[1:5, ] %>% t()
    
data = data.frame(E, G, IND, B, AGE, SEX, nUMI, MT,
                   gPC1 = meta$gPC1, gPC2 = meta$gPC2, gPC3 = meta$gPC3, gPC4 = meta$gPC4, gPC5 = meta$gPC5,
                   expPC1 = expPC[,1], expPC2 = expPC[,2], expPC3 = expPC[,3], expPC4 = expPC[,4], expPC5 = expPC[,5])
data$G = as.numeric(as.character(data$G))

message('Fit Full NBME model')
ptm <- proc.time() # start the stopwatch!

full_model <- lme4::glmer.nb(formula = E ~ G + (1 | IND) + (1 | B) +
                             AGE + SEX + nUMI + MT + 
                             gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
                             expPC1 + expPC2 + expPC3 + expPC4 + expPC5, 
                             data = data, nAGQ = 0, control = glmerControl(optimizer = "nloptwrap"))
out = summary(full_model)$coeff
colnames(out) <- c("Estimate","Std.Error","zvalue","pval")

message('Fit Null NBME model')
null_model <- lme4::glmer.nb(formula = E ~ (1 | IND) + (1 | B) + # no genotype term
                             AGE + SEX + nUMI + MT + 
                             gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
                             expPC1 + expPC2 + expPC3 + expPC4 + expPC5, 
                             data = data, nAGQ = 0, control = glmerControl(optimizer = "nloptwrap"))

message('Run likelihood ratio test')
model_lrt <- lrtest(null_model, full_model)
proc.time() - ptm # end the stopwatch!

message('Save results')
outprefix = '/data/srlab1/jkang/hla2023/eqtl_single_cell/1_sceQTL_nostate/OneK1K_NBME_'
out = data.frame(gene = gene, snp = snp, term = row.names(out), out, 
                 NB_theta = getME(full_model, "glmer.nb.theta"),          # NB overdispersion parameter
                 lrt_ChiSq = model_lrt[2, 4], lrt_pval = model_lrt[2, 5]) # Save LRT p-val for G effect

write.csv(out, paste0(outprefix, paste('result', cell_type, gene, snp, sep = '_'), '.csv'), quote = F)
saveRDS(data, paste0(outprefix, paste('data', cell_type, gene, snp, sep = '_'), '.rds'))