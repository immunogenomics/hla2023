#!/usr/bin/Rscript

# Makes single-cell eQTL plots for OneK1K-defined embedding results, including:
# - UMAP colored estimated beta_total for each cell
# - UMAP colored by quintiles of beta_total
# - Boxplot showing eQTL effect (each dot is 1 sample), by genotype for top and bottom quintile of cells (two versions: free y-axis and shared y-axis)

args = commandArgs(trailingOnly=TRUE)

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
})

fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width, repr.plot.res = 200)
}

## TODO: add individual points
## Add rev for quantile green yellow (quantile_rank --> quantile)
## Add REF/ALT

c = args[1] #'Myeloid'
gene = args[2] #'HLA-DQA1'
variant = args[3] #'rs3104413'
nquantiles = as.numeric(args[4]) #5

# Path to data and output from NBME model
prefix_base = '/data/srlab1/jkang/hla2023/eqtl_single_cell/3_sceQTL_state_interaction_OneK1K/'
prefix_data = paste0(prefix_base, 'OneK1K_NBME_10PCs_data_')
prefix_res = paste0(prefix_base, 'OneK1K_NBME_10PCs_result_')

# Path to reference OneK1K object
prefix_ref = '/data/srlab1/jkang/hla2023/symphony/'

# Path to plotting output folder
prefix_plots = paste0(prefix_base, 'standard_plots/OneK1K')

###############
# Functions
###############

model_quantiles_nbme = function(data, nquantiles) {
    out = NULL
    for (i in 1:nquantiles) {
        datai = data %>% filter(quantile_rank == i)
        full_model <- lme4::glmer.nb(formula = E ~ G + (1 | IND) + (1 | B) +
                             AGE + SEX + nUMI + MT + 
                             gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
                             expPC1 + expPC2 + expPC3 + expPC4 + expPC5, 
                             nAGQ=0, data=datai, control = glmerControl(optimizer = "nloptwrap"))
        outi = summary(full_model)$coeff %>% as.data.frame()
        colnames(outi) <- c("Estimate","Std.Error","zvalue","pval")
        outi$quantile = paste0('Quantile_', i)
        outi$term = rownames(outi)
        out = rbind(out, outi)
    }
    out$quantile = paste(state, as.factor(out$quantile))
    out$quantile = factor(out$quantile)
    outG = out %>% filter(term == 'G')
    outG$pval_format = formatC(outG$pval, format = "e", digits = 2) # scientific notation
    return(outG)
}

plot_boxplot_by_quantile = function(mydata, outG) {
    p = ggplot(mydata, aes(x = G_round, y = mean_E_lognorm)) +
            geom_boxplot(outlier.size = 0.2) + 
            geom_text(data = outG, aes(x = 1, y = 1.1 * max(mydata$mean_E_lognorm), 
                               label = paste('beta[NBME]', "==", round(Estimate, 2))), parse = TRUE) +
            geom_text(data = outG, aes(x = 2.5, y = 1.1 * max(mydata$mean_E_lognorm),
                               label = paste("P", "==", pval_format)), parse = TRUE) +
            theme_classic() + facet_wrap(~ quantile, ncol = length(unique(mydata$quantile))) + 
            xlab(paste0(variant, ' dosage')) + ylab(paste('Mean normalized', gene)) 
    return(p)
}

# Plot only a certain quantile
plot_boxplot_quantile = function(mydata, outG, q) {
    mydata_q = mydata %>% filter(quantile_rank == q)
    outG_q = outG %>% filter(endsWith(as.character(quantile), as.character(q)))
    
    p = ggplot(mydata_q, aes(x = G_round, y = mean_E_lognorm)) +
            geom_boxplot(outlier.size = 0.2) + 
            geom_text(data = outG_q, aes(x = 1, y = 1.1 * max(mydata$mean_E_lognorm), 
                               label = paste('beta[NBME]', "==", round(Estimate, 2))), parse = TRUE) +
            geom_text(data = outG_q, aes(x = 2.5, y = 1.1 * max(mydata$mean_E_lognorm),
                               label = paste("P", "==", pval_format)), parse = TRUE) +
            theme_classic() + 
            xlab(paste0(variant, ' dosage')) + ylab(paste('Mean normalized', gene)) 
    return(p)
}

# Plot only a certain quantile, free y-axis
plot_boxplot_quantile_v2 = function(mydata, outG, q) {
    mydata_q = mydata %>% filter(quantile_rank == q)
    outG_q = outG %>% filter(endsWith(as.character(quantile), as.character(q)))
    
    p = ggplot(mydata_q, aes(x = G_round, y = mean_E_lognorm)) +
            geom_boxplot(outlier.size = 0.2) + 
            geom_text(data = outG_q, aes(x = 1, y = max(1.1 * mydata_q$mean_E_lognorm), 
                               label = paste('beta[NBME]', "==", round(Estimate, 2))), parse = TRUE) +
            geom_text(data = outG_q, aes(x = 2.5, y = max(1.1 * mydata_q$mean_E_lognorm),
                               label = paste("P", "==", pval_format)), parse = TRUE) +
            theme_classic() + 
            xlab(paste0(variant, ' dosage')) + ylab(paste('Mean normalized', gene)) 
    return(p)
}

###############
# Read in data
###############

# Read in data
res = read.csv(paste0(prefix_res, paste(c, gene, variant, sep = '_'), '.csv'))
data = readRDS(paste0(prefix_data, paste(c, gene, variant, sep = '_'), '.rds'))
ref = readRDS(paste0(prefix_ref, 'OneK1K_', c, '_batch2_reference.rds'))

# Find most significant interaction hPC
state = (res %>% filter(startsWith(term, 'G:')) %>% slice(which.min(abs(pval))) %>% dplyr::select(term) %>% as.character() %>% 
        strsplit(split = 'G:') %>% unlist())[2]
print(paste('Most significant cell state interaction is:', state))

# Make the expression by individual and quantile rank
data$quantile_rank = ntile(data[, state], nquantiles)
data$E_lognorm = as.numeric(ref$exp_norm_HLA[gene, ])
data$B = paste0('OneK1K_', ref$meta_data$Batch)

data$G_round = round(data$G) # round to 0, 1, 2
mydata = data %>% group_by(IND, quantile_rank) %>% # take mean of all cells in a quantile per individual
         dplyr::summarise(mean_E_lognorm = mean(E_lognorm),
                          G_round = mean(G_round),
                          quantile = mean(quantile_rank))
mydata$G_round = as.factor(mydata$G_round)
mydata$quantile = paste0(state, ' Quantile_', as.factor(mydata$quantile))
mydata$quantile = factor(mydata$quantile)

###############
# Run NBME model on each quantile, plot boxplots
###############

# Run the NBME model on each quantile separately
outG = model_quantiles_nbme(data, nquantiles)

fig.size(3, 14)
pdf(paste(prefix_plots, c, gene, variant, state, 'Q', nquantiles, '.pdf', sep = '_'), width = 14 * (nquantiles / 5), height = 3)
p = plot_boxplot_by_quantile(mydata, outG)
print(p)
dev.off()

# save a version with factor levels reversed
outG$quantile = forcats::fct_rev(outG$quantile) # optional
mydata$quantile = forcats::fct_rev(mydata$quantile) # optional

pdf(paste(prefix_plots, c, gene, variant, state, 'Q', nquantiles, 'reversed.pdf', sep = '_'), width = 14 * (nquantiles / 5), height = 3)
p = plot_boxplot_by_quantile(mydata, outG)
print(p)
dev.off()

## save ggplot object and data
plottingdata = list()
plottingdata[['mydata']] = mydata
plottingdata[['outG']] = outG

###############
# Calculate and plot single-cell betas on UMAP
###############

# Load eQTL results, CV scores, and UMAP coordinates
sceqtl_results <- res
cv_scores <- data[, c('harmony1', 'harmony2', 'harmony3', 'harmony4', 'harmony5',
                      'harmony6', 'harmony7', 'harmony8', 'harmony9', 'harmony10')]
umap_res <- ref$umap$embedding %>% as.data.frame()

betas = sceqtl_results %>% filter(grepl("^G", term)) %>% dplyr::select(Estimate) %>% unlist
sc_betas = t(data.frame(X1 = rowSums(sweep(cbind(1, cv_scores[,1:10]), MARGIN=2, betas, `*`))))
umap_res$sc_betas = sc_betas %>% as.numeric()
umap_res = cbind(umap_res, t(ref$Z_corr[ 1:10,]))

plottingdata[['umap_res_betas']] = umap_res

fig.size(4, 4.5)
pdf(paste(prefix_plots, c, gene, variant, 'est_betas.pdf', sep = '_'), width = 4.5, height = 4)
p = umap_res %>% 
    sample_frac(1L) %>%
    ggplot() +
        geom_point_rast(aes(x = UMAP1, y = UMAP2, col = sc_betas), size = 0.2) +
        scale_color_gradientn(colours = c("#3c438f", "#7781f2", "lightgrey", "#f2c777", "#d68d06")) +
        theme_void() + labs(col = expression(beta[total]))
print(p)
dev.off()

# Do version with colors reversed
fig.size(4, 4.5)
pdf(paste(prefix_plots, c, gene, variant, 'est_betas_reversedColors.pdf', sep = '_'), width = 4.5, height = 4)
p = umap_res %>% 
    sample_frac(1L) %>%
    ggplot() +
        geom_point_rast(aes(x = UMAP1, y = UMAP2, col = sc_betas), size = 0.2) +
        scale_color_gradientn(colours = rev(c("#3c438f", "#7781f2", "lightgrey", "#f2c777", "#d68d06"))) +
        theme_void() + labs(col = expression(beta[total]))
print(p)
dev.off()

###############
# Plot by bucket of estimated beta-interaction
###############

# Calculate single-cell betas (for interaction terms only!)
intbetas = sceqtl_results %>% filter(grepl("^G:", term)) %>% dplyr::select(Estimate) %>% unlist
sc_intbetas = t(data.frame(X1 = rowSums(sweep(cv_scores[,1:10], MARGIN=2, intbetas, `*`))))

state = 'Est_Int_Beta'
data$sc_intbetas = as.data.frame(t(sc_intbetas))$X1

data$quantile_rank = ntile(data[, 'sc_intbetas'], nquantiles)
mydata = data %>% group_by(IND, quantile_rank) %>% # take mean of all cells in a quantile per individual
         dplyr::summarise(mean_E_lognorm = mean(E_lognorm),
                          G_round = mean(G_round),
                          quantile = mean(quantile_rank))
mydata$G_round = as.factor(mydata$G_round)
mydata$quantile = paste0(state, ' Quantile_', mydata$quantile) %>% factor()

out = model_quantiles_nbme(data, nquantiles)
outG = out %>% filter(term == 'G')
outG$pval_format = formatC(outG$pval, format = "e", digits = 2) # scientific notation
outG$quantile = factor(outG$quantile, levels = paste0(state, ' Quantile_', 1:nquantiles))
mydata$quantile = factor(mydata$quantile, levels = paste0(state, ' Quantile_', 1:nquantiles))

plottingdata[['mydata_estintbeta']] = mydata
plottingdata[['outG_estintbeta']] = outG

pdf(paste(prefix_plots, c, gene, variant, state, 'Q', nquantiles, '.pdf', sep = '_'), width = 14 * (nquantiles / 5), height = 3)
p = ggplot(mydata, aes(x = G_round, y = mean_E_lognorm)) +
    geom_boxplot(outlier.size = 0.2) + 
    # Changed to 1.1 to prevent going off the page
    geom_text(data = outG, aes(x = 1.1, y = max(mydata$mean_E_lognorm), 
                               label = paste('beta[NBME]', "==", round(Estimate, 2))), parse = TRUE, size = 3.5) +
    geom_text(data = outG, aes(x = 2.5, y = max(mydata$mean_E_lognorm),
                               label = paste("P", "==", pval_format)), parse = TRUE, size = 3.5) +
    theme_classic() + facet_wrap(~ quantile, ncol = nquantiles) + 
    xlab(paste0(variant, ' dosage')) + ylab(paste('Mean normalized', gene))
print(p)
dev.off()

## Plot max and min quantiles

# Plot max quantile only
pdf(paste(prefix_plots, c, gene, variant, state, 'Q1', 'of', 
          paste0('Q', nquantiles), 'only.pdf', sep = '_'), , width = 3.5, height = 3)
p = plot_boxplot_quantile_v2(mydata, outG, 1)
print(p)
dev.off()
p

# Plot max quantile only
pdf(paste(prefix_plots, c, gene, variant, state, paste0('Q', nquantiles), 'of', 
          paste0('Q', nquantiles), 'only.pdf', sep = '_'), , width = 3.5, height = 3)
p = plot_boxplot_quantile_v2(mydata, outG, nquantiles)
print(p)
dev.off()
p

# Plot quantiles of intbeta on UMAP
pdf(paste(prefix_plots, c, gene, variant, state, 'Q', nquantiles, '_UMAP.pdf', sep = '_'), width = 4.65, height = 4)
umap_res$quantile_rank = factor(data$quantile_rank)
p = cbind(umap_res, ref$meta_data) %>% 
    sample_frac(1L) %>%
    ggplot() +
        geom_point_rast(aes(x = UMAP1, y = UMAP2, col = quantile_rank), size = 0.2) +
        theme_void() + scale_color_brewer(palette = 'YlGn') +
        guides(colour = guide_legend(override.aes = list(size=4))) +
        labs('Quantile')
print(p)
dev.off()

plottingdata[['umap_res_intbetas']] = umap_res
saveRDS(plottingdata, paste(prefix_plots, 'plottingdata', c, gene, variant, '.rds', sep  = '_')) 