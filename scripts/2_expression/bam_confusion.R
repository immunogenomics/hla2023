#!/usr/bin/Rscript
suppressPackageStartupMessages({
    library(Rsamtools)
    library(tidyr)
    library(ggVennDiagram)
    library(rtracklayer)
    library(Biostrings)
    library(Matrix)
    library(tidyverse)})

fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width, repr.plot.res = 200)
}

args = commandArgs(trailingOnly=TRUE)
sample = args[1] # sample name
results_noPers_dir = args[2] #'/data/srlab2/jkang/scHLA/noPers_final/Randolph2021_NewPanel/STARsolo_results/'
results_pers_dir = args[3] #'/data/srlab2/jkang/scHLA/personalized_final/Randolph2021_NewPanel/STARsolo_results/'

df = as.data.frame(matrix(nrow = 8, ncol = 3))
colnames(df) = c('gene', 'start', 'end')

annotations = rtracklayer::import('/data/srlab1/amber_joyce/filtered_gtf/gencode.v38.annotation.filtered.gtf')
classical = c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DPA1', 'HLA-DPB1')

for (i in 1:length(classical)) {
    classical_annotation = annotations[which(annotations$gene_name == classical[i] &
                                          annotations$type == 'gene'),]
    start = start(ranges(classical_annotation))
    end = end(ranges(classical_annotation))
    df[i, 'gene'] = classical[i]
    df[i, 'start'] = start
    df[i, 'end'] = end
}
rownames(df) = df$gene

message('Reading in bam files')
#########################
# Read in BAM files
#########################
# Read in noPers HLA
param = ScanBamParam(what = c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "seq"))
noPers_bam = scanBam(paste0(results_noPers_dir, sample, '_GeneFull_Ex50pAS/', sample, '_GeneFull_Ex50pAS_HLA.bam'), param=param)
sapply(noPers_bam[[1]], class)
noPers_bam = as.data.frame(noPers_bam) # convert to data frame

# Read in noPers unmapped
noPers_unmapped = scanBam(paste0(results_noPers_dir, sample, '_GeneFull_Ex50pAS/', sample, '_GeneFull_Ex50pAS_unmapped.bam'), param=param)
sapply(noPers_unmapped[[1]], class)
noPers_unmapped = as.data.frame(noPers_unmapped) # convert to data frame
noPers_unmapped_reads = noPers_unmapped$qname
noPers_unmapped_reads = data.frame(qname = noPers_unmapped_reads, gene = 'Unmapped', pipeline = 'noPers')
rm(noPers_unmapped)

# Read in Pers HLA
pers_bam = scanBam(paste0(results_pers_dir, sample, '_GeneFull_Ex50pAS/', sample, '_GeneFull_Ex50pAS_HLA.bam'), param=param)
sapply(pers_bam[[1]], class)
pers_bam = as.data.frame(pers_bam)  # convert to data frame

# Read in pers unmapped
pers_unmapped = scanBam(paste0(results_pers_dir, sample, '_GeneFull_Ex50pAS/', sample, '_GeneFull_Ex50pAS_unmapped.bam'), param=param)
sapply(pers_unmapped[[1]], class)
pers_unmapped = as.data.frame(pers_unmapped) # convert to data frame
pers_unmapped_reads = pers_unmapped$qname
pers_unmapped_reads = data.frame(qname = pers_unmapped_reads, gene = 'Unmapped', pipeline = 'Pers')
rm(pers_unmapped)

# Get list of chromosome names
chr_names = pers_bam$rname %>% unique() %>% as.character()

message('Getting alleles')
# Get alleles
A_alleles = chr_names[which(startsWith(chr_names, 'IMGT_A'))]
B_alleles = chr_names[which(startsWith(chr_names, 'IMGT_B'))]
C_alleles = chr_names[which(startsWith(chr_names, 'IMGT_C'))]
DRB1_alleles = chr_names[which(startsWith(chr_names, 'IMGT_DRB1'))]
DQA1_alleles = chr_names[which(startsWith(chr_names, 'IMGT_DQA1'))]
DQB1_alleles = chr_names[which(startsWith(chr_names, 'IMGT_DQB1'))]
DPA1_alleles = chr_names[which(startsWith(chr_names, 'IMGT_DPA1'))]
DPB1_alleles = chr_names[which(startsWith(chr_names, 'IMGT_DPB1'))]

pers_bam_A = pers_bam %>% subset(rname %in% A_alleles)
pers_bam_B = pers_bam %>% subset(rname %in% B_alleles)
pers_bam_C = pers_bam %>% subset(rname %in% C_alleles)
pers_bam_DQA1 = pers_bam %>% subset(rname %in% DQA1_alleles)
pers_bam_DQB1 = pers_bam %>% subset(rname %in% DQB1_alleles)
pers_bam_DRB1 = pers_bam %>% subset(rname %in% DRB1_alleles)
pers_bam_DPA1 = pers_bam %>% subset(rname %in% DPA1_alleles)
pers_bam_DPB1 = pers_bam %>% subset(rname %in% DPB1_alleles)


message('Getting reads')
# Separate out the reads mapping to each gene in pers
pers_A_reads = unique(pers_bam$qname[which(pers_bam$rname == A_alleles[1])], pers_bam$qname[which(pers_bam$rname == A_alleles[2])])
pers_A_reads = data.frame(qname = pers_A_reads, gene = 'HLA-A', pipeline = 'Pers')

pers_B_reads = unique(pers_bam$qname[which(pers_bam$rname == B_alleles[1])], pers_bam$qname[which(pers_bam$rname == B_alleles[2])])
pers_B_reads = data.frame(qname = pers_B_reads, gene = 'HLA-B', pipeline = 'Pers')

pers_C_reads = unique(pers_bam$qname[which(pers_bam$rname == C_alleles[1])], pers_bam$qname[which(pers_bam$rname == C_alleles[2])])
pers_C_reads = data.frame(qname = pers_C_reads, gene = 'HLA-C', pipeline = 'Pers')

pers_DPA1_reads = unique(pers_bam$qname[which(pers_bam$rname == DPA1_alleles[1])], pers_bam$qname[which(pers_bam$rname == DPA1_alleles[2])])
pers_DPA1_reads = data.frame(qname = pers_DPA1_reads, gene = 'HLA-DPA1', pipeline = 'Pers')

pers_DPB1_reads = unique(pers_bam$qname[which(pers_bam$rname == DPB1_alleles[1])], pers_bam$qname[which(pers_bam$rname == DPB1_alleles[2])])
pers_DPB1_reads = data.frame(qname = pers_DPB1_reads, gene = 'HLA-DPB1', pipeline = 'Pers')

pers_DQA1_reads = unique(pers_bam$qname[which(pers_bam$rname == DQA1_alleles[1])], pers_bam$qname[which(pers_bam$rname == DQA1_alleles[2])])
pers_DQA1_reads = data.frame(qname = pers_DQA1_reads, gene = 'HLA-DQA1', pipeline = 'Pers')

pers_DQB1_reads = unique(pers_bam$qname[which(pers_bam$rname == DQB1_alleles[1])], pers_bam$qname[which(pers_bam$rname == DQB1_alleles[2])])
pers_DQB1_reads = data.frame(qname = pers_DQB1_reads, gene = 'HLA-DQB1', pipeline = 'Pers')

pers_DRB1_reads = unique(pers_bam$qname[which(pers_bam$rname == DRB1_alleles[1])], pers_bam$qname[which(pers_bam$rname == DRB1_alleles[2])])
pers_DRB1_reads = data.frame(qname = pers_DRB1_reads, gene = 'HLA-DRB1', pipeline = 'Pers')

pers_otherMHC_reads = unique(pers_bam$qname[which(pers_bam$rname == 'chr6')])
pers_otherMHC_reads = data.frame(qname = pers_otherMHC_reads, gene = 'Other MHC', pipeline = 'Pers')

# Separate out the reads mapping to each gene in noPers
noPers_A_reads = (noPers_bam %>% subset(pos >= df['HLA-A', 'start'] & pos <= df['HLA-A', 'end']) %>% unique())$qname
noPers_A_reads = data.frame(qname = noPers_A_reads, gene = 'HLA-A', pipeline = 'noPers')

noPers_B_reads = (noPers_bam %>% subset(pos >= df['HLA-B', 'start'] & pos <= df['HLA-B', 'end']) %>% unique())$qname
noPers_B_reads = data.frame(qname = noPers_B_reads, gene = 'HLA-B', pipeline = 'noPers')

noPers_C_reads = (noPers_bam %>% subset(pos >= df['HLA-C', 'start'] & pos <= df['HLA-C', 'end']) %>% unique())$qname
noPers_C_reads = data.frame(qname = noPers_C_reads, gene = 'HLA-C', pipeline = 'noPers')

noPers_DPA1_reads = (noPers_bam %>% subset(pos >= df['HLA-DPA1', 'start'] & pos <= df['HLA-DPA1', 'end']) %>% unique())$qname
noPers_DPA1_reads = data.frame(qname = noPers_DPA1_reads, gene = 'HLA-DPA1', pipeline = 'noPers')

noPers_DPB1_reads = (noPers_bam %>% subset(pos >= df['HLA-DPB1', 'start'] & pos <= df['HLA-DPB1', 'end']) %>% unique())$qname
noPers_DPB1_reads = data.frame(qname = noPers_DPB1_reads, gene = 'HLA-DPB1', pipeline = 'noPers')

noPers_DQA1_reads = (noPers_bam %>% subset(pos >= df['HLA-DQA1', 'start'] & pos <= df['HLA-DQA1', 'end']) %>% unique())$qname
noPers_DQA1_reads = data.frame(qname = noPers_DQA1_reads, gene = 'HLA-DQA1', pipeline = 'noPers')

noPers_DQB1_reads = (noPers_bam %>% subset(pos >= df['HLA-DQB1', 'start'] & pos <= df['HLA-DQB1', 'end']) %>% unique())$qname
noPers_DQB1_reads = data.frame(qname = noPers_DQB1_reads, gene = 'HLA-DQB1', pipeline = 'noPers')

noPers_DRB1_reads = (noPers_bam %>% subset(pos >= df['HLA-DRB1', 'start'] & pos <= df['HLA-DRB1', 'end']) %>% unique())$qname
noPers_DRB1_reads = data.frame(qname = noPers_DRB1_reads, gene = 'HLA-DRB1', pipeline = 'noPers')

noPers_otherMHC_reads = noPers_bam$qname[which(!(noPers_bam$pos >= df['HLA-A', 'start'] & noPers_bam$pos <= df['HLA-A', 'end']) &
                                                      !(noPers_bam$pos >= df['HLA-B', 'start'] & noPers_bam$pos <= df['HLA-B', 'end']) &
                                                      !(noPers_bam$pos >= df['HLA-C', 'start'] & noPers_bam$pos <= df['HLA-C', 'end']) &
                                                      !(noPers_bam$pos >= df['HLA-DPA1', 'start'] & noPers_bam$pos <= df['HLA-DPA1', 'end']) &
                                                      !(noPers_bam$pos >= df['HLA-DPB1', 'start'] & noPers_bam$pos <= df['HLA-DPB1', 'end']) &
                                                      !(noPers_bam$pos >= df['HLA-DQA1', 'start'] & noPers_bam$pos <= df['HLA-DQA1', 'end']) &
                                                      !(noPers_bam$pos >= df['HLA-DQB1', 'start'] & noPers_bam$pos <= df['HLA-DQB1', 'end']) &
                                                      !(noPers_bam$pos >= df['HLA-DRB1', 'start'] & noPers_bam$pos <= df['HLA-DRB1', 'end'] ))] %>% unique()
noPers_otherMHC_reads = data.frame(qname = noPers_otherMHC_reads, gene = 'Other MHC', pipeline = 'noPers')

pers_reads = rbind(pers_A_reads, pers_B_reads, pers_C_reads, 
                   pers_DPA1_reads, pers_DPB1_reads, pers_DQA1_reads, pers_DQB1_reads, pers_DRB1_reads, 
                   pers_otherMHC_reads, pers_unmapped_reads) %>% mutate(Sample = sample)
noPers_reads = rbind(noPers_A_reads, noPers_B_reads, noPers_C_reads, 
                     noPers_DPA1_reads, noPers_DPB1_reads, noPers_DQA1_reads, noPers_DQB1_reads, noPers_DRB1_reads, 
                     noPers_otherMHC_reads, noPers_unmapped_reads) %>% mutate(Sample = sample)

message('Making results')
results = NULL
for (gene in c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRB1')) {
    reads_to_follow = pers_reads$qname[which(pers_reads$gene == gene)] # read names mapping to pers gene
    noPers_reads_fromgene = noPers_reads %>% filter(qname %in% reads_to_follow)
    t = table(noPers_reads_fromgene$gene) %>% as.data.frame() %>% mutate(pers_gene = gene)
    results = plyr::rbind.fill(results, t)
}
colnames(results) = c('noPers_read_alignment', 'n', 'pers_read_alignment')
results$noPers_read_alignment = factor(results$noPers_read_alignment, levels = 
                                       c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRB1', 'Other MHC', 'Unmapped'))

saveRDS(results, paste0('/data/srlab1/jkang/hla/notebooks/220906_bam_confusion/', sample, '_noPers_read_alignments.rds'))

results = NULL
for (gene in c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRB1')) {
    reads_to_follow = noPers_reads$qname[which(noPers_reads$gene == gene)] # read names mapping to noPers gene
    pers_reads_fromgene = pers_reads %>% filter(qname %in% reads_to_follow)
    t = table(pers_reads_fromgene$gene) %>% as.data.frame() %>% mutate(noPers_gene = gene)
    results = plyr::rbind.fill(results, t)
}
colnames(results) = c('pers_read_alignment', 'n', 'noPers_read_alignment')
results$pers_read_alignment = factor(results$pers_read_alignment, levels = 
                                       c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRB1', 'Other MHC', 'Unmapped'))

saveRDS(results, paste0('/data/srlab1/jkang/hla/notebooks/220906_bam_confusion/', sample, '_pers_read_alignments.rds'))
