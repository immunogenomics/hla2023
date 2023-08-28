#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

# Perform sc-eQTL interaction analysis in each dataset separately,
# Using the hPCs defined by mapping blood states onto tissue reference.

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

cell_type = args[1] #'T'
gene = args[2]      # 'HLA-DQA1'
snp = args[3]       #'rs3104371'

# Read in tissue hPCs for state interaction test
refquery = readRDS(paste0('/data/srlab1/jkang/hla2023/symphony/RefQuery_', cell_type,'_mapBloodOntoTissue.rds'))
colnames(refquery$Z_corr) = refquery$meta_data$Cell

############
#  OneK1K  #
#######################################################################
message('Running sc-eQTL model on OneK1K dataset')
ref_path = paste0('/data/srlab1/jkang/hla2023/symphony/OneK1K_', cell_type, '_batch2_reference.rds') # path to ref object
ref = readRDS(ref_path)
# Load genotype dosages
geno = readRDS('/data/srlab1/jkang/hla2023/data/sampleXdosage/OneK1K_sampleXdosage_final.rds')

# Load raw UMI counts matrix and metadata and gene expression PCs from Symphony object
meta = ref$meta_data
pca_res = ref$Z_orig # pre-Harmony PCs
harmonypca_res = refquery$Z_corr[, meta$Cell] # Harmony PCs comes from refquery object

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
harmonyPC = harmonypca_res[1:10, ] %>% t()
    
data = data.frame(E, G, IND, B, AGE, SEX, nUMI, MT,
                   gPC1 = meta$gPC1, gPC2 = meta$gPC2, gPC3 = meta$gPC3, gPC4 = meta$gPC4, gPC5 = meta$gPC5, dataset = dataset,
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
                             data = data, nAGQ = 0, control = glmerControl(optimizer = "nloptwrap"))
out = summary(full_model)$coeff
colnames(out) <- c("Estimate","Std.Error","zvalue","pval")

message('Fit Null NBME model')
null_model <- lme4::glmer.nb(formula = E ~ G + (1 | IND) + (1 | B) +
                             AGE + SEX + nUMI + MT + 
                             gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
                             expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + 
                             harmony1 + harmony2 + harmony3 + harmony4 + harmony5 +
                             harmony6 + harmony7 + harmony8 + harmony9 + harmony10, # no G*state interaction terms
                             data = data, nAGQ = 0, control = glmerControl(optimizer = "nloptwrap"))

message('Run likelihood ratio test')
model_lrt <- lrtest(null_model, full_model)

message('Save results')
outprefix = '/data/srlab1/jkang/hla2023/eqtl_single_cell/2_sceQTL_state_interaction_tissue_hPCs/OneK1K_mapToTissue_NBME_10PCs_'
out = data.frame(gene = gene, snp = snp, dataset = 'OneK1K',  # mapped to tissue
                 cell_type = cell_type, term = row.names(out), out, 
                 lrt_pval = model_lrt[2, 5], lrt_ChiSq = model_lrt[2, 4]) # Save LRT p-val for G*state effect
write.csv(out, paste0(outprefix, paste('result', cell_type, gene, snp, sep = '_'), '.csv'), quote = F)
saveRDS(data, paste0(outprefix, paste('data', cell_type, gene, snp, sep = '_'), '.rds'))

#################
#  Randolph2021 #
#######################################################################
message('Running sc-eQTL model on Randolph2021 dataset')
ref_path = paste0('/data/srlab1/jkang/hla2023/symphony/Randolph_', cell_type, '_batch2_reference.rds') # path to ref object
ref = readRDS(ref_path)
# Load genotype dosages
geno = readRDS('/data/srlab1/jkang/hla2023/data/sampleXdosage/Randolph2021_sampleXdosage_final.rds')

# Load raw UMI counts matrix and metadata and gene expression PCs from Symphony object
meta = ref$meta_data
pca_res = ref$Z_orig # pre-Harmony PCs
harmonypca_res = refquery$Z_corr[, meta$Cell] # Harmony PCs comes from refquery object

message('Make dataframe')
# Make data frame of variables for model
E = as.numeric(ref$exp_HLA[gene, ]) %>% round() # round to get integers
G = t(geno[, snp])[, as.character(meta$Geno_ID)] # changed to Geno_ID for Randolph nomenclature
IND = as.factor(meta$Geno_ID)
B = paste0('Randolph_', meta$Batch)
AGE = scale(meta$Age)
SEX = meta$Female
nUMI = scale(log(meta$nUMI))
MT = scale(meta$percent.mito)
dataset = meta$dataset
expPC = pca_res[1:5, ] %>% t()
harmonyPC = harmonypca_res[1:10, ] %>% t()
    
data = data.frame(E, G, IND, B, AGE, SEX, nUMI, MT,
                   gPC1 = meta$gPC1, gPC2 = meta$gPC2, gPC3 = meta$gPC3, gPC4 = meta$gPC4, gPC5 = meta$gPC5, dataset = dataset,
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
                             data = data, nAGQ = 0, control = glmerControl(optimizer = "nloptwrap"))
out = summary(full_model)$coeff
colnames(out) <- c("Estimate","Std.Error","zvalue","pval")

message('Fit Null NBME model')
null_model <- lme4::glmer.nb(formula = E ~ G + (1 | IND) + (1 | B) +
                             AGE + SEX + nUMI + MT + 
                             gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
                             expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + 
                             harmony1 + harmony2 + harmony3 + harmony4 + harmony5 +
                             harmony6 + harmony7 + harmony8 + harmony9 + harmony10, # no G*state interaction terms
                             data = data, nAGQ = 0, control = glmerControl(optimizer = "nloptwrap"))

message('Run likelihood ratio test')
model_lrt <- lrtest(null_model, full_model)

message('Save results')
outprefix = '/data/srlab1/jkang/hla2023/eqtl_single_cell/2_sceQTL_state_interaction_tissue_hPCs/Randolph2021_mapToTissue_NBME_10PCs_'
out = data.frame(gene = gene, snp = snp, dataset = 'Randolph2021',  # mapped to tissue
                 cell_type = cell_type, term = row.names(out), out, 
                 lrt_pval = model_lrt[2, 5], lrt_ChiSq = model_lrt[2, 4]) # Save LRT p-val for G*state effect
write.csv(out, paste0(outprefix, paste('result', cell_type, gene, snp, sep = '_'), '.csv'), quote = F)
saveRDS(data, paste0(outprefix, paste('data', cell_type, gene, snp, sep = '_'), '.rds'))


############
#  AMP2RA  #
#######################################################################
message('Running sc-eQTL model on AMP2RA dataset')
ref_path = paste0('/data/srlab1/jkang/hla2023/symphony/AMP2RA_', cell_type, '_batch2_reference.rds') # path to ref object
ref = readRDS(ref_path)
# Load genotype dosages
geno = readRDS('/data/srlab1/jkang/hla2023/data/sampleXdosage/AMP2RA_sampleXdosage_final.rds')

# Load raw UMI counts matrix and metadata and gene expression PCs from Symphony object
meta = ref$meta_data
pca_res = ref$Z_orig # pre-Harmony PCs
harmonypca_res = refquery$Z_corr[, meta$Cell] # Harmony PCs comes from refquery object 

message('Make dataframe')
# Make data frame of variables for model
E = as.numeric(ref$exp_HLA[gene, ]) %>% round() # round to get integers
G = t(geno[, snp])[, as.character(meta$Sample)]
IND = as.factor(meta$Geno_ID)
#No Batch for AMP2RA
AGE = scale(meta$Age)
SEX = meta$Female
nUMI = scale(log(meta$nUMI))
MT = scale(meta$percent.mito)
dataset = meta$dataset
expPC = pca_res[1:5, ] %>% t()
harmonyPC = harmonypca_res[1:10, ] %>% t()
    
data = data.frame(E, G, IND, AGE, SEX, nUMI, MT,
                   gPC1 = meta$gPC1, gPC2 = meta$gPC2, gPC3 = meta$gPC3, gPC4 = meta$gPC4, gPC5 = meta$gPC5, dataset = dataset,
                   expPC1 = expPC[,1], expPC2 = expPC[,2], expPC3 = expPC[,3], expPC4 = expPC[,4], expPC5 = expPC[,5], 
                   harmony1 = harmonyPC[, 1], harmony2 = harmonyPC[, 2], harmony3 = harmonyPC[, 3], 
                   harmony4 = harmonyPC[, 4], harmony5 = harmonyPC[, 5], 
                   harmony6 = harmonyPC[, 6], harmony7 = harmonyPC[, 7], harmony8 = harmonyPC[, 8], 
                   harmony9 = harmonyPC[, 9], harmony10 = harmonyPC[, 10])
data$G = as.numeric(as.character(data$G))

message('Fit Full NBME model')
full_model <- lme4::glmer.nb(formula = E ~ G + (1 | IND) + #(1 | B) +
                             AGE + SEX + nUMI + MT + 
                             gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
                             expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + 
                             harmony1 + G:harmony1 + harmony2 + G:harmony2 + harmony3 + G:harmony3 + 
                             harmony4 + G:harmony4 + harmony5 + G:harmony5 + 
                             harmony6 + G:harmony6 + harmony7 + G:harmony7 + harmony8 + G:harmony8 + 
                             harmony9 + G:harmony9 + harmony10 + G:harmony10, 
                             data = data, nAGQ = 0, control = glmerControl(optimizer = "nloptwrap"))
out = summary(full_model)$coeff
colnames(out) <- c("Estimate","Std.Error","zvalue","pval")

message('Fit Null NBME model')
null_model <- lme4::glmer.nb(formula = E ~ G + (1 | IND) + #(1 | B) +
                             AGE + SEX + nUMI + MT + 
                             gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
                             expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + 
                             harmony1 + harmony2 + harmony3 + harmony4 + harmony5 +
                             harmony6 + harmony7 + harmony8 + harmony9 + harmony10, # no G*state interaction terms
                             data = data, nAGQ = 0, control = glmerControl(optimizer = "nloptwrap"))

message('Run likelihood ratio test')
model_lrt <- lrtest(null_model, full_model)

message('Save results')
outprefix = '/data/srlab1/jkang/hla2023/eqtl_single_cell/2_sceQTL_state_interaction_tissue_hPCs/AMP2RA_mapToTissue_NBME_10PCs_'
out = data.frame(gene = gene, snp = snp, dataset = 'AMP2RA',  # tissue
                 cell_type = cell_type, term = row.names(out), out, 
                 lrt_pval = model_lrt[2, 5], lrt_ChiSq = model_lrt[2, 4]) # Save LRT p-val for G*state effect
write.csv(out, paste0(outprefix, paste('result', cell_type, gene, snp, sep = '_'), '.csv'), quote = F)
saveRDS(data, paste0(outprefix, paste('data', cell_type, gene, snp, sep = '_'), '.rds'))

message('All done!')