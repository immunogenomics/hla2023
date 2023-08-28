#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

## Calculates power to detect genotype effects across a range of allele frequencies and effect sizes

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
    library(hexbin)
    library(ggpubr)
    library(patchwork)
})

source('/data/srlab1/jkang/hla2023/scripts/utils.R')

# simulation parameters
p = as.numeric(args[1]) #0.05 # allele frequency
genotype_effect = as.numeric(args[2]) #0.8

coef_results = c()
sig_results = c()

prefix = '/data/srlab1/jkang/hla2023/eqtl_single_cell/'
simulated_counts = readRDS(paste0(prefix, '5_sceQTL_simulations/simulated_counts.rds'))
data = readRDS(paste0(prefix, '5_sceQTL_simulations/data_predictors.rds'))

set.seed(1)
for (i in 1:length(simulated_counts)) {
    print(i)
    ### Simulate the genotype dosages at desired allele frequency

    # H-W equilibrium:
    # p = 0.5 (frequency of allele A)
    # q = 1-p (frequence of allele a)
    # p^2 + 2*pq + q^2 = 1 # H-W equilibrium
    n = length(unique(data$IND)) # num individuals

    dosage_2 = p^2
    dosage_1 = 2*p*(1-p)
    dosage_0 = (1-p)^2
    genos_sim = c(rep(2, n*dosage_2), rep(1, n*dosage_1), rep(0, n*dosage_0))
    extra = n - length(genos_sim)
    genos_sim = c(genos_sim, sample(c(0, 1, 2), extra))

    ### "Spike in" prespecified genotype effect to simulated expression (E_sim)

    # Add the simulated genotypes to the data for modeling
    sim_G = data.frame(cbind(unique(as.character(data$IND)), sample(genos_sim)))
    colnames(sim_G) = c('IND', 'G_sim')
    data_sim = left_join(data, sim_G)
    data_sim$G_sim = as.numeric(data_sim$G_sim)

    # Generate the simulated expression data with "spiked in" genotype effect
    data_sim$E_sim = exp(log(simulated_counts[[i]]) + genotype_effect * data_sim$G_sim) %>% round()

    ### Fit full model on simulated data

    full_model_sim <- lme4::glmer.nb(formula = E_sim ~ G_sim + (1 | IND) + (1 | B) +
                             AGE + SEX + nUMI + MT + 
                             gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
                             expPC1 + expPC2 + expPC3 + expPC4 + expPC5, 
                             data = data_sim, nAGQ = 0, control = glmerControl(optimizer = "nloptwrap"))
    null_model_sim <- lme4::glmer.nb(formula = E_sim ~ (1 | IND) + (1 | B) +
                             AGE + SEX + nUMI + MT + 
                             gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
                             expPC1 + expPC2 + expPC3 + expPC4 + expPC5, 
                             data = data_sim, nAGQ = 0, control = glmerControl(optimizer = "nloptwrap"))

    model_lrt <- lrtest(null_model_sim, full_model_sim)
    coef_results[i] = summary(full_model_sim)$coeff['G_sim', 'Estimate']
    sig_results[i] = model_lrt[2, 5]
}

### Save results
res = cbind(coef_results, sig_results) %>% as.data.frame()
res$AF = p
res$Gbeta = genotype_effect

message('Save results')
saveRDS(res, paste0(prefix, '5_sceQTL_simulations/', 'res_freq_', 
                    p, '_Gbeta_', genotype_effect, '.rds'))