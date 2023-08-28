#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

# Makes single-cell eQTL plots for tissue embedding results for OneK1K and AMP2RA datasets, including:
# - UMAP colored estimated beta_total for each cell
# - UMAP colored by quintiles of beta_total
# - Boxplot showing eQTL effect (each dot is 1 sample), by genotype for top and bottom quintile of cells (two versions: free y-axis and shared y-axis)

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
    library(ggplot2)
    library(ggrastr)
    library(RColorBrewer)
    library(forcats)
    library(ggrepel)
    library(symphony)
    library(irlba)
    library(patchwork)
})
source('/data/srlab1/jkang/hla2023/scripts/utils.R')

cell_type = args[1] #'T'
gene = args[2] # 'HLA-DQA1'
variant = args[3] #'rs3104371'
nquantiles = as.numeric(args[4]) #5

sc_betas_combined = NULL
plots = list()
myorder = c(0, 0.25, 0.5, 0.75, 1)

### Useful functions
# Plot eQTL boxplot for a certain quantile - mean(log2(UMI+1))
plot_boxplot_quantile_v3 = function(mydata, outG, q) {
    mydata_q = mydata %>% filter(quantile_rank == q)
    outG_q = outG %>% filter(endsWith(as.character(quantile), as.character(q)))
    
    p = ggplot(mydata_q, aes(x = G_round, y = mean_E)) +
            geom_boxplot(outlier.size = 0.2) + 
            geom_jitter_rast(color="black", size=0.3, alpha=0.9, width = 0.2) +
            theme_classic(base_size = 13) + 
            xlab(variant) + ylab(paste('Mean log2(UMI+1)', gene)) + expand_limits(y = 0)
    return(p)
}

# Plot eQTL boxplot for a certain quantile - mean(log2(UMI+1)) - except keeps the same Y-axis across quantiles
plot_boxplot_quantile_v4 = function(mydata, outG, q) {
    Ymax = max(mydata$mean_E) + 0.001
    
    mydata_q = mydata %>% filter(quantile_rank == q)
    outG_q = outG %>% filter(endsWith(as.character(quantile), as.character(q)))
    
    p = ggplot(mydata_q, aes(x = G_round, y = mean_E)) +
            geom_boxplot(outlier.size = 0.2) + 
            geom_text(data = outG_q, aes(x = 1, y = max(1.1 * mydata$mean_E), 
                               label = paste('beta[NBME]', "==", round(Estimate, 2))), parse = TRUE, size = 4.5) +
            geom_text(data = outG_q, aes(x = 2.5, y = max(1.1 * mydata$mean_E),
                               label = paste("P", "==", pval_format)), parse = TRUE, size = 4.5) +
            geom_jitter_rast(color="black", size=0.3, alpha=0.9, width = 0.2) +
            theme_classic(base_size = 13) + 
            xlab(variant) + ylab(paste('Mean log2(UMI+1)', gene)) + 
            ylim(0, 1.1 * Ymax)
    return(p)
}

# Run NBME model on each quantile separately
model_quantiles_nbme = function(data, nquantiles) {
    out = NULL
    for (i in 1:nquantiles) {
        datai = data %>% filter(quantile_rank == i)
        full_model <- lme4::glmer.nb(formula = E ~ G + (1 | IND) + (1 | B) +
                             AGE + SEX + nUMI + MT + 
                             gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
                             expPC1 + expPC2 + expPC3 + expPC4 + expPC5, 
                             nAGQ=0, data=datai, control = glmerControl(optimizer = "nloptwrap"))
        
        null_model <- lme4::glmer.nb(formula = E ~ (1 | IND) + (1 | B) +
                             AGE + SEX + nUMI + MT + 
                             gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
                             expPC1 + expPC2 + expPC3 + expPC4 + expPC5, 
                             nAGQ=0, data=datai, control = glmerControl(optimizer = "nloptwrap"))
        
        lrt_pval = lrtest(full_model, null_model)[2, 5]
        outi = summary(full_model)$coeff %>% as.data.frame()
        colnames(outi) <- c("Estimate","Std.Error","zvalue","pval")
        outi$quantile = paste0('Quantile_', i)
        outi$term = rownames(outi)
        outi$lrt_pval = lrt_pval
        out = rbind(out, outi)
    }
    out$quantile = paste(state, as.factor(out$quantile))
    out$quantile = factor(out$quantile)
    outG = out %>% filter(term == 'G')
    outG$pval_format = formatC(outG$lrt_pval, format = "e", digits = 2) # scientific notation
    return(outG)
}

# Run NBME model on each quantile separately (for AMP2RA)
model_quantiles_nbme_noBatch = function(data, nquantiles) {
    out = NULL
    for (i in 1:nquantiles) {
        datai = data %>% filter(quantile_rank == i)
        full_model <- lme4::glmer.nb(formula = E ~ G + (1 | IND) + # (1 | B) +
                             AGE + SEX + nUMI + MT + 
                             gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
                             expPC1 + expPC2 + expPC3 + expPC4 + expPC5, 
                             nAGQ=0, data=datai, control = glmerControl(optimizer = "nloptwrap"))
        
        null_model <- lme4::glmer.nb(formula = E ~ (1 | IND) + # (1 | B) +
                             AGE + SEX + nUMI + MT + 
                             gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
                             expPC1 + expPC2 + expPC3 + expPC4 + expPC5, 
                             nAGQ=0, data=datai, control = glmerControl(optimizer = "nloptwrap"))
        lrt_pval = lrtest(full_model, null_model)[2, 5]
        
        outi = summary(full_model)$coeff %>% as.data.frame()
        colnames(outi) <- c("Estimate","Std.Error","zvalue","pval")
        outi$quantile = paste0('Quantile_', i)
        outi$term = rownames(outi)
        outi$lrt_pval = lrt_pval
        out = rbind(out, outi)
    }
    out$quantile = paste(state, as.factor(out$quantile))
    out$quantile = factor(out$quantile)
    outG = out %>% filter(term == 'G')
    outG$pval_format = formatC(outG$lrt_pval, format = "e", digits = 2) # scientific notation
    return(outG)
}

message('Starting...')
# Path to reference object
prefix_ref = '/data/srlab1/jkang/hla2023/symphony/'

# Path for plot and plotting data
prefix_plot = '/data/srlab1/jkang/hla2023/eqtl_single_cell/2_sceQTL_state_interaction_tissue_hPCs/standard_plots/'

tissueRef_query = readRDS(paste0('/data/srlab1/jkang/hla2023/symphony/RefQuery_', cell_type,'_mapBloodOntoTissue.rds'))

for (dataset in c('AMP2RA', 'OneK1K')) {
    message('Calc sc betas for ', dataset)
    
    # Path to data and output from sc-eQTL model
    prefix_data = paste0('/data/srlab1/jkang/hla2023/eqtl_single_cell/2_sceQTL_state_interaction_tissue_hPCs/', 
                         dataset, '_mapToTissue_NBME_10PCs_data_')
    prefix_res = paste0('/data/srlab1/jkang/hla2023/eqtl_single_cell/2_sceQTL_state_interaction_tissue_hPCs/', 
                        dataset, '_mapToTissue_NBME_10PCs_result_')

    # Read in data
    res = read.csv(paste0(prefix_res, paste(cell_type, gene, variant, sep = '_'), '.csv'))
    data = readRDS(paste0(prefix_data, paste(cell_type, gene, variant, sep = '_'), '.rds'))
    
    dataset_ref = readRDS(paste0(prefix_ref, dataset, '_', cell_type, '_batch2_reference.rds'))

    # Fix rownames for indexing
    rownames(tissueRef_query$umap$embedding) = tissueRef_query$meta_data$Cell
    colnames(tissueRef_query$Z_corr) = tissueRef_query$meta_data$Cell
    
    # Load eQTL results, CV scores, and UMAP coordinates
    sceqtl_results <- res
    cv_scores <- data[, paste0('harmony', 1:10)]
    
    umap_res <- tissueRef_query$umap$embedding[dataset_ref$meta_data$Cell, ] %>% as.data.frame()

    # Calculate single-cell betas
    betas = sceqtl_results %>% filter(grepl("^G", term)) %>% dplyr::select(Estimate) %>% unlist
    sc_betas = t(data.frame(X1 = rowSums(sweep(cbind(1, cv_scores[,1:10]), MARGIN=2, betas, `*`))))
    umap_res$sc_betas = sc_betas %>% as.numeric()
    umap_res = cbind(umap_res, t(tissueRef_query$Z_corr[1:10, dataset_ref$meta_data$Cell]))
    umap_res$dataset = dataset
    sc_betas_combined = rbind(sc_betas_combined, umap_res)
    
    p = umap_res %>% 
        sample_frac(1L) %>%
        ggplot(aes(x = UMAP1, y = UMAP2, z = sc_betas)) +
            theme_void(base_size = 14) + labs(col = expression(beta[total])) +
            geom_point_rast(aes(col = sc_betas), size = 0.2) +
            scale_color_gradientn(colours = c("#3c438f", "#7781f2", "lightgrey", "#f2c777", "#d68d06"), 
                              values = myorder)+labs(col = expression(beta[total]))
    plots[[dataset]] = p
    pdf(paste(prefix_plot, cell_type, gene, variant, dataset, 'est_betas.pdf', sep = '_'), width = 4.5, height = 3.8)
    print(p)
    dev.off()
    
    # Plot a version with reversed colors
    pdf(paste(prefix_plot, cell_type, gene, variant, dataset, 
              'est_betas_revColors.pdf', sep = '_'), width = 4.5, height = 3.8)
    p = p + scale_color_gradientn(colours = c("#3c438f", "#7781f2", "lightgrey", "#f2c777", "#d68d06"), 
                              values = rev(myorder)) +labs(col = expression(beta[total]))
    print(p)
    dev.off()
    plots[[paste0(dataset, '_revColors')]] = p
    
    ### Calculate single-cell betas (for interaction terms only!)
    data$E = as.numeric(dataset_ref$exp_HLA[gene, ])
    data$E_lognorm = as.numeric(dataset_ref$exp_norm_HLA[gene, ])
    data$G_round = round(data$G)
    intbetas = sceqtl_results %>% filter(grepl("^G:", term)) %>% dplyr::select(Estimate) %>% unlist
    sc_intbetas = t(data.frame(X1 = rowSums(sweep(cv_scores[,1:10], MARGIN=2, intbetas, `*`))))

    state = 'Est_Int_Beta'
    data$sc_intbetas = as.data.frame(t(sc_intbetas))$X1

    ### Run the NBME model on each quantile of interaction beta
    data$quantile_rank = ntile(data[, 'sc_intbetas'], nquantiles)
    mydata = data %>% group_by(IND, quantile_rank) %>% # take mean of all cells in a quantile per individual
         dplyr::summarise(mean_E_lognorm = mean(E_lognorm),
                          mean_E = mean(log2(E + 1)),
                          G_round = mean(G_round),
                          quantile = mean(quantile_rank))
    mydata$G_round = as.factor(mydata$G_round)
    mydata$quantile = paste0(state, ' Quantile_', mydata$quantile) %>% factor()

    if (dataset != 'AMP2RA') {
        outG = model_quantiles_nbme(data, nquantiles)
    } else {
        outG = model_quantiles_nbme_noBatch(data, nquantiles)
    }
    outG$quantile = factor(outG$quantile, levels = paste0(state, ' Quantile_', 1:nquantiles))
    mydata$quantile = factor(mydata$quantile, levels = paste0(state, ' Quantile_', 1:nquantiles))
    
    ## Change REF/ALT to nucleotides
    geno_df = readRDS('/data/srlab1/jkang/hla2023/data/sampleXdosage/four_cohorts_variant_info_final.rds')
    geno_df$POS = as.numeric(geno_df$POS)
    
    myREF = geno_df[which(geno_df$ID == variant), 'REF']
    myALT = geno_df[which(geno_df$ID == variant), 'ALT']
    mydata$G_round = plyr::mapvalues(mydata$G_round, c('0', '1', '2'), 
                                      c(paste0(myREF, '/', myREF), paste0(myREF, '/', myALT), 
                                        paste0(myALT, '/', myALT)), warn_missing = TRUE)
    
    ## Plot eQTL boxplots
    pdf(paste(prefix_plot, cell_type, gene, variant, dataset, state, 'Q', nquantiles, '_log2UMI.pdf', sep = '_'), 
        width = 14 * (nquantiles / 5), height = 3)
    p = ggplot(mydata, aes(x = G_round, y = mean_E)) +
        geom_boxplot(outlier.size = 0.2) + 
        # Changed to 1.1 to prevent going off the page
        geom_text(data = outG, aes(x = 1.1, y = max(mydata$mean_E), 
                               label = paste('beta[NBME]', "==", round(Estimate, 2))), parse = TRUE, size = 3.5) +
        geom_text(data = outG, aes(x = 2.5, y = max(mydata$mean_E),
                               label = paste("P", "==", pval_format)), parse = TRUE, size = 3.5) +
        geom_jitter(color="black", size=0.3, alpha=0.9, width = 0.2) +
        theme_classic(base_size = 14) + facet_wrap(~ quantile, ncol = nquantiles) + 
        xlab(variant) + ylab(paste('Mean log2(UMI+1)', gene))
    print(p)
    dev.off()

    ## Plot max and min quantiles with different axes
    # Plot min quantile only
    pdf(paste(prefix_plot, cell_type, gene, variant, dataset, state, 'Q1', 'of', 
          paste0('Q', nquantiles), 'only_log2UMI.pdf', sep = '_'), width = 3.5, height = 3)
    p = plot_boxplot_quantile_v3(mydata, outG, 1)
    print(p)
    dev.off()
    p

    # Plot max quantile only
    pdf(paste(prefix_plot, cell_type, gene, variant, dataset, state, paste0('Q', nquantiles), 'of', 
          paste0('Q', nquantiles), 'only_log2UMI.pdf', sep = '_'), width = 3.5, height = 3)
    p = plot_boxplot_quantile_v3(mydata, outG, nquantiles)
    print(p)
    dev.off()
    p
    
    ## Plot max and min quantiles with shared axes
    
    # Plot min quantile only
    pdf(paste(prefix_plot, cell_type, gene, variant, dataset, state, 'Q1', 'of', 
          paste0('Q', nquantiles), 'only_log2UMI_sharedYax.pdf', sep = '_'), width = 3.5, height = 3)
    p = plot_boxplot_quantile_v4(mydata, outG, 1)
    print(p)
    dev.off()
    p

    # Plot max quantile only
    pdf(paste(prefix_plot, cell_type, gene, variant, dataset, state, paste0('Q', nquantiles), 'of', 
          paste0('Q', nquantiles), 'only_log2UMI_sharedYax.pdf', sep = '_'), width = 3.5, height = 3)
    p = plot_boxplot_quantile_v4(mydata, outG, nquantiles)
    print(p)
    dev.off()
    p

    ## Plot quantiles of interaction beta on UMAP (various shades of green)
    pdf(paste(prefix_plot, cell_type, gene, variant, dataset, state, 'Q', nquantiles, 
              '_UMAP.pdf', sep = '_'), width = 4.65, height = 3.8)
    umap_res$quantile_rank = factor(data$quantile_rank)
    p = umap_res %>% 
        sample_frac(1L) %>%
        ggplot() +
            geom_point_rast(aes(x = UMAP1, y = UMAP2, col = quantile_rank), size = 0.2) +
            theme_void(base_size = 14) + scale_color_brewer(palette = 'YlGn') +
            guides(colour = guide_legend(override.aes = list(size=4))) +
            labs(col = expression(atop(beta[total], quantile)))
    print(p)
    dev.off()
    plots[[paste0(dataset, '_quantiles')]] = p
    
    pdf(paste(prefix_plot, cell_type, gene, variant, dataset, state, 'Q', nquantiles, 
              '_UMAP_revColors.pdf', sep = '_'), width = 4.65, height = 3.8)
    umap_res$quantile_rank = factor(data$quantile_rank)
    p = umap_res %>% 
        sample_frac(1L) %>%
        ggplot() +
            geom_point_rast(aes(x = UMAP1, y = UMAP2, col = quantile_rank), size = 0.2) +
            theme_void(base_size = 14) + scale_color_brewer(palette = 'YlGn', direction = -1) +
            guides(colour = guide_legend(override.aes = list(size=4))) +
            labs(col = expression(atop(beta[total], quantile)))
    print(p)
    dev.off()
    plots[[paste0(dataset, '_quantiles_revColors')]] = p
}

min_UMAP1 = min(sc_betas_combined[which(sc_betas_combined$dataset %in% c('AMP2RA', 'OneK1K')), 'UMAP1'])
max_UMAP1 = max(sc_betas_combined[which(sc_betas_combined$dataset %in% c('AMP2RA', 'OneK1K')), 'UMAP1'])
min_UMAP2 = min(sc_betas_combined[which(sc_betas_combined$dataset %in% c('AMP2RA', 'OneK1K')), 'UMAP2'])
max_UMAP2 = max(sc_betas_combined[which(sc_betas_combined$dataset %in% c('AMP2RA', 'OneK1K')), 'UMAP2'])

pdf(paste(prefix_plot, cell_type, gene, variant, 'AMP2RA_OneK1K_est_betas.pdf', sep = '_'), width = 9.5, height = 3.8)
q1 = (plots[['AMP2RA']] + ylim(min_UMAP2, max_UMAP2) + xlim(min_UMAP1, max_UMAP1) ) + plot_spacer() +
     (plots[['OneK1K']] + ylim(min_UMAP2, max_UMAP2) + xlim(min_UMAP1, max_UMAP1) ) + 
     plot_layout(widths = c(4.5, 0.5, 4.5))
print(q1)
dev.off()

pdf(paste(prefix_plot, cell_type, gene, variant, 'AMP2RA_OneK1K_est_betas_revColors.pdf', sep = '_'), width = 9.5, height = 3.8)
q1 = (plots[['AMP2RA_revColors']] + ylim(min_UMAP2, max_UMAP2) + xlim(min_UMAP1, max_UMAP1) ) + plot_spacer() +
     (plots[['OneK1K_revColors']] + ylim(min_UMAP2, max_UMAP2) + xlim(min_UMAP1, max_UMAP1) ) + 
     plot_layout(widths = c(4.5, 0.5, 4.5))
print(q1)
dev.off()

pdf(paste(prefix_plot, cell_type, gene, variant, 'AMP2RA_OneK1K_Qs.pdf', sep = '_'), width = 9.5, height = 3.8)
q1 = (plots[['AMP2RA_quantiles']] + ylim(min_UMAP2, max_UMAP2) + xlim(min_UMAP1, max_UMAP1) ) + plot_spacer() +
     (plots[['OneK1K_quantiles']] + ylim(min_UMAP2, max_UMAP2) + xlim(min_UMAP1, max_UMAP1) ) + 
     plot_layout(widths = c(4.65, 0.2, 4.65))
print(q1)
dev.off()

pdf(paste(prefix_plot, cell_type, gene, variant, 'AMP2RA_OneK1K_Qs_revColors.pdf', sep = '_'), width = 9.5, height = 3.8)
q1 = (plots[['AMP2RA_quantiles_revColors']] + ylim(min_UMAP2, max_UMAP2) + xlim(min_UMAP1, max_UMAP1) ) + plot_spacer() +
     (plots[['OneK1K_quantiles_revColors']] + ylim(min_UMAP2, max_UMAP2) + xlim(min_UMAP1, max_UMAP1) ) + 
     plot_layout(widths = c(4.65, 0.2, 4.65))
print(q1)
dev.off()

saveRDS(sc_betas_combined, paste0(paste(prefix_plot, 'data', cell_type, gene, variant, 'est_betas.rds', sep = '_')))