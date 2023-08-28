#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

# Tests for calibration of cell state interaction analysis using NBME model, for a given cell type, gene, and variant, using the OneK1K dataset and OneK1K-defined hPC embedding.
# Significance is defined as LRT between full model and reduced model without cell state
# Permutes cell state across donors randomly, runs nperms permutations. 

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

set.seed(0)

message('Read in single-cell data')
cell_type = args[1] #'T'
ref_path = paste0('/data/srlab1/jkang/hla2023/symphony/OneK1K_', cell_type, '_batch2_reference.rds') # path to ref object
gene = args[2] # 'HLA-DQA1'
snp = args[3] #'rs3104371'
ref = readRDS(ref_path)
nperms = 1000

# Load raw UMI counts matrix and metadata and gene expression PCs from Symphony object
meta = ref$meta_data
pca_res = ref$Z_orig # pre-Harmony PCs
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

results = NULL
waldPs = NULL
for (i in 1:nperms) {
    
    harmonyPC = harmonypca_res[1:10, ] %>% t()
    
    # Shuffle the harmony PC block among all cells (across donors)
    harmonyPC = harmonyPC[sample(rownames(harmonyPC)), 1:10]
    
    data = data.frame(E, G, B, IND, AGE, SEX, nUMI, MT,
                gPC1 = meta$gPC1, gPC2 = meta$gPC2, gPC3 = meta$gPC3, gPC4 = meta$gPC4, gPC5 = meta$gPC5,
                expPC1 = expPC[,1], expPC2 = expPC[,2], expPC3 = expPC[,3], expPC4 = expPC[,4], expPC5 = expPC[,5], 
                harmony1 = harmonyPC[, 1], harmony2 = harmonyPC[, 2], harmony3 = harmonyPC[, 3], 
                harmony4 = harmonyPC[, 4], harmony5 = harmonyPC[, 5], 
                harmony6 = harmonyPC[, 6], harmony7 = harmonyPC[, 7], harmony8 = harmonyPC[, 8], 
                harmony9 = harmonyPC[, 9], harmony10 = harmonyPC[, 10])
    data$G = as.numeric(as.character(data$G))

    message('Fit Full NBME model')
    full_model <- lme4::glmer.nb(formula = E ~ G + (1 | IND) + (1 | B) +
                        AGE + SEX + nUMI + MT + 
                        gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
                        expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + 
                        harmony1 + G:harmony1 + harmony2 + G:harmony2 + harmony3 + G:harmony3 + 
                        harmony4 + G:harmony4 + harmony5 + G:harmony5 + 
                        harmony6 + G:harmony6 + harmony7 + G:harmony7 + harmony8 + G:harmony8 + 
                        harmony9 + G:harmony9 + harmony10 + G:harmony10, 
                        nAGQ=0, data= data, control = glmerControl(optimizer = "nloptwrap"))
    out = summary(full_model)$coeff
    colnames(out) <- c("Estimate","Std.Error","zvalue","pval")

    message('Fit Null NBME model')
    null_model <- lme4::glmer.nb(formula = E ~ G + (1 | IND) + (1 | B) +
                        AGE + SEX + nUMI + MT + 
                        gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
                        expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + 
                        harmony1 + harmony2 + harmony3 + harmony4 + harmony5 +
                        harmony6 + harmony7 + harmony8 + harmony9 + harmony10, # no G*state interaction terms
                        nAGQ=0, data= data, control = glmerControl(optimizer = "nloptwrap"))

    message('Run likelihood ratio test')
    model_lrt <- lrtest(null_model, full_model)
    
    result = data.frame(gene = gene, snp = snp, perm_i = i, lrt_ChiSq = model_lrt[2, 4], lrt_pval = model_lrt[2, 5])
    results = rbind(result, results)
    out = data.frame(out)
    out$perm_i = i
    waldPs = rbind(waldPs, out)
}

message('Save results')
outprefix = '/data/srlab1/jkang/hla2023/eqtl_single_cell/4_sceQTL_permutation_OneK1K/'
saveRDS(results, paste0(outprefix, paste('NBME_OneK1K_hPCS_', cell_type, gene, snp, 
                                         nperms, 'perms_acrossDonors', sep = '_'), '.rds'))
saveRDS(waldPs, paste0(outprefix, paste('NBME_OneK1K_hPCS_', cell_type, gene, snp, 
                                        nperms, 'WaldPs_acrossDonors', sep = '_'), '.rds'))