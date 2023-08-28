#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

# This script will estimate the percent of variance of each gene's expression can be explained by cell state (hPCs).

# Load libraries
suppressPackageStartupMessages({
    library(tidyr)
    library(dplyr)
    library(plyr)
    library(Matrix)
    library(Matrix.utils)
    library(singlecellmethods)
    library(irlba)
    library(symphony)
    library(ggplot2)
    library(corrplot)
    library(pheatmap)
    library(viridis)
    library(ggrastr)
    library(RColorBrewer)
    library(dplyr)
    library(purrr)
    library(ggplot2)
    library(cowplot)
    library(ggpubr)
    library(lme4)
    library(ggrepel)
    require(MuMIn)
})

# Inputs and constants
dataset = args[1]
prefix = '/data/srlab1/jkang/hla2023/'
dataset_res = NULL
for (c in c('B_plasma', 'Myeloid', 'T')) {
    
    meta = readRDS(paste0(prefix, 'eqtl_pseudobulk/0_data/cell_meta_', dataset, '_', c ,'.rds'))
    meta$X = NULL
    sample_meta = read.csv(paste0(prefix, 'data/meta/sample_meta_', dataset, '_completeHLA.csv'))
    
    meta = left_join(meta, sample_meta)
    exp = readRDS(paste0(prefix, 'eqtl_pseudobulk/0_data/exp_pers_', dataset, '_', c,'.rds'))
    
    # Get the hPCs
    refquery = readRDS(paste0(prefix, 'symphony/RefQuery_', c,'_mapBloodOntoTissue.rds'))
    dataset_name = dataset
    if (dataset == 'Smillie') dataset_name = 'Smillie2019'
    if (dataset == 'Randolph') dataset_name = 'Randolph2021'
    Z_symphony = refquery$Z_corr[, which(refquery$meta$dataset == dataset_name)]
    
    ref = readRDS(paste0(prefix, 'symphony/', dataset, '_', c, '_batch2_reference.rds'))
    meta = cbind(meta, t(Z_symphony), t(ref$Z_orig))
    
    if (dataset %in% c('AMP2RA', 'Smillie')) { # no batch random effect
        for (gene in c('A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1')) {
            # Add gene counts
            name = paste0('HLA-', gene)
            meta = cbind(meta, exp[name, ])
            colnames(meta)[length(colnames(meta))] = paste0('HLA-', gene)
            full <- lme4::glmer.nb(formula = round(get(name)) ~ (1 | Sample) + Age + Female + 
                    scale(log(nUMI)) + scale(percent.mito) + 
                    harmony_1 + harmony_2 + harmony_3 + harmony_4 + harmony_5 +
                            harmony_6 + harmony_7 + harmony_8 + harmony_9 + harmony_10 +
                    gPC1 + gPC2 + gPC3 + gPC4 + gPC5, nAGQ=0, 
                    data = meta, control = glmerControl(optimizer = "nloptwrap"))
            full_rsq = r.squaredGLMM(full)['delta', 'R2m']
        
            nostate <- lme4::glmer.nb(formula = round(get(name)) ~ (1 | Sample) + Age + Female + 
                    scale(log(nUMI)) + scale(percent.mito) + 
                    gPC1 + gPC2 + gPC3 + gPC4 + gPC5, nAGQ=0, 
                    data = meta, control = glmerControl(optimizer = "nloptwrap"))
            nostate_rsq = r.squaredGLMM(nostate)['delta', 'R2m']
        
            cellstate_rsq = full_rsq - nostate_rsq
            out = data.frame(cell_type = c, gene = gene, full_rsq = full_rsq, nostate_rsq = nostate_rsq,
                         cellstate_rsq = cellstate_rsq)
            dataset_res = rbind(dataset_res, out)   
        }
    } else {
        for (gene in c('A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1')) {
            # Add gene counts
            name = paste0('HLA-', gene)
            meta = cbind(meta, exp[name, ])
            colnames(meta)[length(colnames(meta))] = paste0('HLA-', gene)
            full <- lme4::glmer.nb(formula = round(get(name)) ~ (1 | Sample) + (1 | Batch) + Age + Female + 
                    scale(log(nUMI)) + scale(percent.mito) + 
                    harmony_1 + harmony_2 + harmony_3 + harmony_4 + harmony_5 +
                            harmony_6 + harmony_7 + harmony_8 + harmony_9 + harmony_10 +
                    gPC1 + gPC2 + gPC3 + gPC4 + gPC5, nAGQ=0, 
                    data = meta, control = glmerControl(optimizer = "nloptwrap"))
            full_rsq = r.squaredGLMM(full)['delta', 'R2m']
        
            nostate <- lme4::glmer.nb(formula = round(get(name)) ~ (1 | Sample) + (1 | Batch) + Age + Female + 
                    scale(log(nUMI)) + scale(percent.mito) + 
                    gPC1 + gPC2 + gPC3 + gPC4 + gPC5, nAGQ=0, 
                    data = meta, control = glmerControl(optimizer = "nloptwrap"))
            nostate_rsq = r.squaredGLMM(nostate)['delta', 'R2m']
        
            cellstate_rsq = full_rsq - nostate_rsq
            out = data.frame(cell_type = c, gene = gene, full_rsq = full_rsq, nostate_rsq = nostate_rsq,
                         cellstate_rsq = cellstate_rsq)
            dataset_res = rbind(dataset_res, out)   
        }
    }
    saveRDS(dataset_res, paste0(prefix, 'fig4/', dataset, '_perc_var_explained.rds'))
}