# Useful functions and variables, including color definitions and gene sets
# Last updated June 2023

fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width, repr.plot.res = 200)
}

mito.genes.coding = c('MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8','MT-ATP6','MT-CO3',
                      'MT-ND3','MT-ND4L','MT-ND4','MT-ND5','MT-ND6','MT-CYB')

# Color definitions for figures
cohort_colors = c("Synovium" = "#92d66b",
                  "Intestine" = "#fac878",
                  "PBMC-cultured" = "#9177CF",
                  "PBMC-blood" = "#84a9e8")
cohort_colors_darker = colorspace::darken(cohort_colors, 0.3)
names(cohort_colors_darker) = names(cohort_colors)

gene_colors = c('HLA-A' = '#f56342',
                'HLA-B' = '#eba83d',
                'HLA-C' = '#d1be15',
                'HLA-DPA1' = '#32ad6e',
                'HLA-DPB1' = '#0f8d96',
                'HLA-DQA1' = '#418ef2',
                'HLA-DQB1' = '#a188db',
                'HLA-DRB1' = '#c76392')

# Function to process SNP2HLA output from new reference
make_imputation_results = function(vcf_dosage_file) {
    vcf = read.vcfR(vcf_dosage_file, verbose = FALSE)
    idx_HLA = grepl("HLA", vcf@fix[, "ID"], fixed = TRUE) # get HLA alleles
    sum(idx_HLA)
    vcf_hla = vcf[idx_HLA]
    vcf_fix_df = as.data.frame(vcf_hla@fix)
    
    # Format columns
    vcf_fix_df = separate(data = vcf_fix_df, col = INFO, into = c("AR2", "DR2", "AF", "IMP"), sep = ";")
    vcf_fix_df$GENE = sub("\\*.*", "", vcf_fix_df$ID)
    vcf_fix_df$AF = sub(".*=", "", vcf_fix_df$AF) %>% as.numeric()
    vcf_fix_df$AR2 = sub(".*=", "", vcf_fix_df$AR2) %>% as.numeric()
    vcf_fix_df$DR2 = sub(".*=", "", vcf_fix_df$DR2) %>% as.numeric()
    vcf_fix_df$ncolons = str_count(vcf_fix_df$ID, ":")
    vcf_fix_df$ndigit = 2 * (vcf_fix_df$ncolons + 1)
    length(unique(vcf_fix_df$ID))
    
    # Make quality plots and tables
    p1 = vcf_fix_df %>% ggplot() + geom_point(aes(x = AF, y = AR2, 
        col = GENE), size = 0.3) + geom_vline(xintercept = 0.01) + 
        ggtitle("All alleles imputation quality (AR2)") + theme_classic() + 
        theme(plot.title = element_text(hjust = 0.5))
    p2 = vcf_fix_df %>% ggplot() + geom_point(aes(x = AF, y = DR2, 
        col = GENE), size = 0.3) + geom_vline(xintercept = 0.01) + 
        ggtitle("All alleles imputation quality (DR2)") + theme_classic() + 
        theme(plot.title = element_text(hjust = 0.5))
    r1 = vcf_fix_df %>% group_by(GENE) %>% 
        summarise(mean = mean(AR2))
    r2 = vcf_fix_df %>% group_by(GENE) %>% 
        summarise(mean = mean(DR2))
    return(list(vcf = vcf, vcf_hla = vcf_hla, vcf_fix_df = vcf_fix_df, 
        p1 = p1, r1 = r1, p2 = p2, r2 = r2))
}

old_new_dataset_names = c("AMP2RA"="Synovium", 
                     "Smillie"="Intestine",
                     "Randolph_NI"="PBMC-cultured",
                     "OneK1K" = "PBMC-blood",
                     "Smillie2019"="Intestine",
                     "Randolph2021"="PBMC-cultured",
                     "Randolph"="PBMC-cultured")

old_new_B_plasma = c("AMP2RA_B-2: IgM+IgD+TCL1A+ naive"="Naive", 
                     "AMP2RA_B-1: CD24++CD27+IgM+ unswitched memory"="Memory",
                     "AMP2RA_B-0: CD24+CD27+CD11b+ switched memory"="Memory",
                     "AMP2RA_B-8: IgG1+IgG3+ plasma" = "Plasma",
                     "AMP2RA_B-7: HLA-DR+IgG+ plasmablast" = "Plasmablast",
                     "AMP2RA_B-3: IgM+IgD+CD1c+ MZ-like" = "MZ-like",
                     "AMP2RA_B-5: CD11c+LAMP1+ ABC" = "ABC",
                     "AMP2RA_B-6: IgM+ plasma" = "Plasma",
                     "AMP2RA_B-4: AICDA+BCL6+ GC-like" = "GC-like",
                     "OneK1K_B naive" = "Naive",
                     "OneK1K_B memory" = "Memory",
                     "OneK1K_B intermediate" = "Intermediate",
                     "OneK1K_Plasmablast" = "Plasmablast")

B_plasma_colors = c("ABC" = "#de1a1a",
                  "GC-like" = "#a80874",
                  "Intermediate" = "#c3b59f",
                  "Memory" = "#247ba0",
                  "MZ-like" = "#d6a2ad",
                  "Naive" = "#b2dbbf", 
                  "Plasma" = "#f2d398", 
                  "Plasmablast" = "#d78521")

old_new_Myeloid =  c('AMP2RA_M-10: DC2' = 'DC2', 
                     'AMP2RA_M-0: MERTK+ SELENOP+ LYVE1+' = "Tissue Macrophage",
                     'AMP2RA_M-7: IL1B+ FCN1+'="Monocyte (Infiltrating)",
                     'AMP2RA_M-9: DC3' = "DC3",
                     'AMP2RA_M-5: C1QA+' = "Monocyte (Infiltrating)",
                     'AMP2RA_M-1: MERTK+ SELENOP+ LYVE1-' = "Tissue Macrophage",
                     'AMP2RA_M-2: MERTK+ S100A8+' = "Tissue Macrophage",
                     'AMP2RA_M-6: STAT1+ CXCL10+' = "Monocyte (Infiltrating)",
                     'AMP2RA_M-4: SPP1+' = "Mono/Mac",
                     'AMP2RA_M-13: pDC' = "pDC",
                     'AMP2RA_M-14: LAMP3+' = "DC (other)",
                     "AMP2RA_M-8: PLCG2+" = "Mono/Mac",
                     "AMP2RA_M-3: MERTK+ HBEGF+" = "Mono/Mac",
                     "AMP2RA_M-12: DC1" = "DC1",
                     'AMP2RA_M-11: DC4' = 'DC4',
                     'OneK1K_CD16 Mono' = 'Monocyte (Blood)',
                     'OneK1K_CD14 Mono' = 'Monocyte (Blood)',
                     'OneK1K_cDC2' = 'DC2',
                     'OneK1K_pDC' = 'pDC',
                     'OneK1K_ASDC' = 'DC (other)',
                     'OneK1K_cDC1' = 'DC1')

Myeloid_colors = c("Monocyte (Blood)" = "#ffe548",
                   "Monocyte (Infiltrating)" = "red",
                   "Mono/Mac" = "maroon",
                   "Tissue Macrophage" = "darkorange",
                   "DC1" = "skyblue",
                   "DC2" = "#3993dd",
                   "DC3" = "#4d967d",
                   "DC4" = "#7dd66f", 
                   "pDC" = "#11356C",
                   "DC (other)" = "tan")


old_new_T =  c('AMP2RA_T-6: CD4+ memory' = 'CD4+ Memory', 
               'AMP2RA_T-5: CD4+ GZMK+ memory' = 'CD4+ Memory',
               'AMP2RA_T-8: CD4+ CD25-high Treg' = 'Treg',
               'AMP2RA_T-22: Vdelta1' = 'gdT',
               'AMP2RA_T-14: CD8+ GZMK+ memory' = 'CD8+ Memory',
               'AMP2RA_T-2: CD4+ IL7R+CCR5+ memory' = 'CD4+ Memory',
               'AMP2RA_T-12: CD4+ GNLY+' = 'CD4+ Cytotoxic',
               'AMP2RA_T-4: CD4+ naive' = 'CD4+ Naive',
               'AMP2RA_T-0: CD4+ IL7R+ memory' = 'CD4+ Memory',
               'AMP2RA_T-3: CD4+ Tfh/Tph' = 'Tph/Tfh',
               'AMP2RA_T-1: CD4+ CD161+ memory' = 'CD4+ Memory',
               'AMP2RA_T-16: CD8+ CD45ROlow/naive' = 'CD8+ Naive',
               'AMP2RA_T-18: Proliferating' = 'Proliferating',
               'AMP2RA_T-20: CD38+' = 'CD4+ Memory',
               'AMP2RA_T-13: CD8+ GZMK/B+ memory' = 'CD8+ Cytotoxic',
               'AMP2RA_T-19: MT-high (low quality)' = 'MT-high',
               'AMP2RA_T-10: CD4+ OX40+NR3C1+' = 'CD4+ Memory',
               'AMP2RA_T-9: CD4+ CD25-low Treg' = 'Treg',
               'AMP2RA_T-15: CD8+ GZMB+/TEMRA' = 'CD8+ Cytotoxic',
               'AMP2RA_T-11: CD4+ CD146+ memory' = 'CD4+ Memory',
               'AMP2RA_T-7: CD4+ Tph' = 'Tph/Tfh',
               'AMP2RA_T-23: Vdelta2' = 'gdT',
               'AMP2RA_T-17: CD8+ activated/NK-like' = 'CD8+ Cytotoxic',
               'OneK1K_CD4 TCM' = 'CD4+ Memory',
               'OneK1K_CD4 Naive' = 'CD4+ Naive',
               'OneK1K_CD8 TEM' = 'CD8+ Cytotoxic',
               'OneK1K_CD8 TCM' = 'CD8+ Memory',
               'OneK1K_CD4 TEM' = 'CD4+ Cytotoxic',
               'OneK1K_Treg' = 'Treg',
               'OneK1K_gdT' = 'gdT',
               'OneK1K_CD8 Naive' = 'CD8+ Naive',
               'OneK1K_dnT' = 'dnT',
               'OneK1K_CD4 CTL' = 'CD4+ Cytotoxic',
               'OneK1K_CD4 Proliferating' = 'Proliferating',
               'OneK1K_CD8 Proliferating' = 'Proliferating')

T_colors = c("CD8+ Memory" = "#eb2a4b",
             "CD4+ Memory" = "#ed95a4",
             "CD8+ Cytotoxic" = "#436bb5",
             "CD4+ Cytotoxic" = "#a6c2f7",
             "CD8+ Naive" = "#2ea361",
             "CD4+ Naive" = "#abedc8",
             "Tph/Tfh" = "#eb7fe9",
             "Treg" = "#e0cd3a",
             "MT-high" = "#54585e",
             "Proliferating" = "#910d3e",
             "gdT" = "#845ad1",
             "dnT" = "#c6e667")

# IFN gene signatures from MSigDB
IFN_alpha_sig = c("ADAR","BATF2","BST2","C1S","CASP1","CASP8","CCRL2","CD47","CD74","CMPK2","CNP",
                "CSF1","CXCL10","CXCL11","DDX60","DHX58","EIF2AK2","ELF1","EPSTI1","FAM125A","FAM46A",
                "FTSJD2","GBP2","GBP4","GMPR","HERC6","IFI27","IFI30","IFI35","IFI44","IFI44L", # rm "HLA-C", "B2M"
                "IFIH1","IFIT2","IFIT3","IFITM1","IFITM2","IFITM3","IL15","IL4R","IL7","IRF1","IRF2",
                "IRF7","IRF9","ISG15","ISG20","LAMP3","LAP3","LGALS3BP","LPAR6","LY6E","MOV10","MX1",
                "NCOA7","NMI","NUB1","OAS1","OASL","OGFR","PARP12","PARP14","PARP9","PLSCR1","PNPT1",
                "PRIC285","PROCR","PSMA3","PSMB8","PSMB9","PSME1","PSME2","RIPK2","RNF31","RSAD2","RTP4",
                "SAMD9","SAMD9L","SELL","SLC25A28","SP110","STAT2","TAP1","TDRD7","TMEM140","TRAFD1","TRIM14",
                "TRIM21","TRIM25","TRIM26","TRIM5","TXNIP","UBA7","UBE2L6","USP18","WARS")

IFN_gamma_sig = c("ADAR","APOL6","ARID5B","ARL4A","AUTS2","B2M","BANK1","BATF2","BPGM","BST2","BTG1",
                "C1R","C1S","CASP1","CASP3","CASP4","CASP7","CASP8","CCL2","CCL5","CCL7","CD274","CD38",
                "CD40","CD69","CD74","CD86","CDKN1A","CFB","CFH","CIITA","CMKLR1","CMPK2","CSF2RB","CXCL10",
                "CXCL11","CXCL9","RIGI","DDX60","DHX58","EIF2AK2","EIF4E3","EPSTI1","FAS","FCGR1A","FGL2","FPR1",
                "CMTR1","GBP4","GBP6","GCH1","GPR18","GZMA","HERC6","HIF1A","HLA-DMA", #"HLA-DQA1","HLA-A","HLA-B", "HLA-DRB1"
                "HLA-G","ICAM1","IDO1","IFI27","IFI30","IFI35","IFI44","IFI44L","IFIH1","IFIT1","IFIT2",
                "IFIT3","IFITM2","IFITM3","IFNAR2","IL10RA","IL15","IL15RA","IL18BP","IL2RB","IL4R","IL6","IL7",
                "IRF1","IRF2","IRF4","IRF5","IRF7","IRF8","IRF9","ISG15","ISG20","ISOC1","ITGB7","JAK2","KLRK1",
                "LAP3","LATS2","LCP2","LGALS3BP","LY6E","LYSMD2","MARCHF1","METTL7B","MT2A","MTHFD2","MVP","MX1",
                "MX2","MYD88","NAMPT","NCOA3","NFKB1","NFKBIA","NLRC5","NMI","NOD1","NUP93","OAS2","OAS3","OASL",
                "OGFR","P2RY14","PARP12","PARP14","PDE4B","PELI1","PFKP","PIM1","PLA2G4A","PLSCR1","PML","PNP","PNPT1",
                "HELZ2","PSMA2","PSMA3","PSMB10","PSMB2","PSMB8","PSMB9","PSME1","PSME2","PTGS2","PTPN1","PTPN2","PTPN6",
                "RAPGEF6","RBCK1","RIPK1","RIPK2","RNF213","RNF31","RSAD2","RTP4","SAMD9L","SAMHD1","SECTM1","SELP",
                "SERPING1","SLAMF7","SLC25A28","SOCS1","SOCS3","SOD2","SP110","SPPL2A","SRI","SSPN","ST3GAL5","ST8SIA4",
                "STAT1","STAT2","STAT3","STAT4","TAP1","TAPBP","TDRD7","TNFAIP2","TNFAIP3","TNFAIP6","TNFSF10","TOR1B",
                "TRAFD1","TRIM14","TRIM21","TRIM25","TRIM26","TXNIP","UBE2L6","UPP1","USP18","VAMP5","VAMP8","VCAM1","WARS1",
                "XAF1","XCL1","ZBP1","ZNFX1")

IFN_sig = union(IFN_alpha_sig, IFN_gamma_sig)

IFN_sig_Davenport = c('HERC5', 'IFI27', 'IRF7', 'ISG15', 'LY6E', 'MX1', 'OAS2', 'OAS3', 'RSAD2', 'USP18', 'GBP5')