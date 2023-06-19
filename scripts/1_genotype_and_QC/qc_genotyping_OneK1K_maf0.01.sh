# @Author: Joyce Kang
# @Date: 2022-03-24
# Adapted from original script by Saori Sakaue & Michelle Curtis
# Genotyping QC for the OneK1K cohort (hg19)

# Genotyping data from Jose Alquicira-Hernandez
DATADIR=/data/srlab/external-data/1K1K_raw/geno/hg19
F="10-final"
FOUT="1K1K"
OUTDIR=/data/srlab/external-data/1K1K_raw/geno/hg19/qc
SCRIPTS=/data/srlab1/jkang/hla/schla/scripts/1_genotype

# Check if its needed to convert to forward strand using snpflip: /data/srlab/external-data/1K1K_raw/geno/hg19/snpflip_output
conda activate env_r3.6
cd /data/srlab/external-data/1K1K_raw/geno/hg19
snpflip --fasta-genome /data/srlab/ssakaue/data/hg19/human_g1k_v37.for_flip.fasta --bim-file $DATADIR/$F.bim

##### Flip the 247,054 SNPs on the reverse strand
conda activate hla_new
module load plink/1.90
plink --bfile $DATADIR/$F --flip snpflip_output/$F.reverse --make-bed  --out $OUTDIR/$FOUT.fw  # Get reverse
#### 496794 variants loaded from .bim file.
#### 975 people (414 males, 561 females) loaded from .fam.
#### --flip: 247054 SNPs flipped, 1 SNP ID not present.

####################
### Deduplication ##
####################

cd /data/srlab1/jkang/hla/schla/scripts
plink --bfile $OUTDIR/$FOUT.fw --missing --out $OUTDIR/$FOUT.fw
python 1_genotype/get_duprem_var.py $OUTDIR/$FOUT.fw
plink --bfile $OUTDIR/$FOUT.fw --exclude $OUTDIR/$FOUT.fw.remdup.snp --make-bed --out $OUTDIR/$FOUT.dedup
#### Removed 682 SNP duplicates
#### 496112 variants and 975 people pass filters and QC

## Remove multiallelic variants
awk '{print $1,$2,$4}' $OUTDIR/$FOUT.dedup.bim|sort|uniq -d > $OUTDIR/dup_pos_rsid.txt
awk '{print $2}' $OUTDIR/dup_pos_rsid.txt > $OUTDIR/dup_pos_rsid.rsids.txt # get list of rsids to exclude
# no multiallelic variants

############
# START QC #
############

# (1) Identification of individuals with discordant sex information --> not needed (skip)
#plink --vcf $DATADIR/$F1 --check-sex --out $OUTPUTDIR/dg_broad
#grep PROBLEM $OUTPUTDIR/dg_broad.sexcheck > $OUTPUTDIR/dg_broad.sexprobs

# (2) Initial SNP QC (prelim for sample QC)
plink --bfile $OUTDIR/$FOUT.dedup --geno 0.1 --hwe 1e-10 --make-bed --out $OUTDIR/$FOUT.1stSNP_QC
#### 0 variants removed due to missing genotype data
#### 0 variants removed due to Hardy-Weinberg exact test 

#############
# Sample QC #
#############

# (3) Identification of individuals with elevated missing data rates
plink --bfile $OUTDIR/$FOUT.1stSNP_QC --missing --out $OUTDIR/$FOUT.1stSNP_QC --allow-no-sex

cd $OUTDIR
# In R
name = '1K1K'
dat = read.table(paste0(name, ".1stSNP_QC.imiss"), header = T)
png(paste0(name, ".1stSNP_QC.imiss.png"))
hist(dat$F_MISS, breaks = 200)
abline(v = 0.01, col = 'red')
dev.off()

plink -bfile $OUTDIR/$FOUT.1stSNP_QC --mind 0.01 --make-bed --out $OUTDIR/$FOUT.1stSNP_QC.mis
#### 3 people removed due to missing genotype data (--mind).
#### IDs written to /data/srlab/external-data/1K1K_raw/geno/hg19/qc/1K1K.1stSNP_QC.mis.irem
######### 505_506, 553_554, 565_566
#### 496112 variants and 972 people pass filters and QC.

# (4) Identification of individuals with outlying heterozygosity rate
plink --bfile $OUTDIR/$FOUT.1stSNP_QC.mis --het --out $OUTDIR/$FOUT.1stSNP_QC.mis --allow-no-sex
### report written to /data/srlab/external-data/1K1K_raw/geno/hg19/qc/1K1K.1stSNP_QC.mis.het

# In R
name = '1K1K'
dat = read.table(paste0(name, ".1stSNP_QC.mis.het") , header=T)
dat$het = (dat$N.NM.-dat$O.HOM.)/dat$N.NM.

# Calculate mean and standard dev
m = mean(dat$het)
print(m) # 0.213
sd = sd(dat$het)
print(sd) # 0.00139
dat$sample = paste(dat$FID, dat$IID, sep = '_')

png(paste0(name, ".1stSNP_QC.mis.het", "_abline.png"))
hist(dat$het, breaks = 200, xlab = "Heterozygosity rate")
abline(v= c(m + (3 * sd),  m - (3 * sd)), col="blue", lwd=2, lty=2, xpd = T) 
# line at +/- 3 x sd
dev.off()
    
# record high het samples
high = dat[dat$het > m + (3 * sd) | dat$het < m - (3 * sd), ]
print(high$sample) # "210_211" "358_359" "512_513" "823_824" "952_953" "975_976"
write.table(high, paste0(name, ".het_outside_3sd.sample.txt"), quote=F, col.names=T, row.names=F)
write.table(high[,1:2], paste0(name, ".het_outside_3sd.sample.fam"), quote=F, col.names=T, row.names=F)

# Do not remove high het
#####################

# (5) Sample relatedness
cd $OUTDIR
HLA_START=24000000
HLA_STOP=36000000

gawk '$1==6 && $4 > '${HLA_START}' && $4 < '${HLA_STOP}'{print $2}' $FOUT.1stSNP_QC.mis.bim | sort > $FOUT.hla.snps
plink --bfile $FOUT.1stSNP_QC.mis --exclude $FOUT.hla.snps --maf 0.05 --indep 50 5 2 --out $FOUT.1stSNP_QC.mis.noHLA.MAF
plink --bfile $FOUT.1stSNP_QC.mis --extract $FOUT.1stSNP_QC.mis.noHLA.MAF.prune.in --genome --out $FOUT.1stSNP_QC.mis.noHLA.MAF
awk '{if($10>0.9)print}' $FOUT.1stSNP_QC.mis.noHLA.MAF.genome > $FOUT.1stSNP_QC.mis.ident.txt
python $SCRIPTS/get_remID.py $FOUT.1stSNP_QC.imiss $FOUT.1stSNP_QC.mis.ident.txt $FOUT.1stSNP_QC.mis.ident.remid.fam
plink --bfile $FOUT.1stSNP_QC.mis --remove $FOUT.1stSNP_QC.mis.ident.remid.fam --make-bed --out $FOUT.1stSNP.1stSamp_QC
cat $FOUT.1stSNP_QC.mis.ident.remid.fam | wc -l # number removed
# 0 individuals removed

##############
# Variant QC #
##############

# Variant QC
plink --bfile $FOUT.1stSNP.1stSamp_QC --missing --out $FOUT.1stSNP.1stSamp_QC

# In R
name = '1K1K'
dat = read.table(paste0(name, ".1stSNP.1stSamp_QC.lmiss"), header=T)
png(paste0(name, ".1stSNP.1stSamp_QC.lmiss.png"))
hist(dat$F_MISS, breaks=200)
abline(v=0.02, col="red")
dev.off()

# Remove missing genotypes
# File names may be incorrect
plink --bfile $FOUT.1stSNP.1stSamp_QC --geno 0.02 --make-bed --out $FOUT.1stSNP.1stSamp_QC.mis
##### 91 variants removed due to missing genotype data (--geno).
##### 496021 variants and 972 people pass filters and QC.


# QC on allele frequency differences with 1KG
LANG=C
bash ${SCRIPTS}/CheckAlleleFreq.sh $FOUT.1stSNP.1stSamp_QC.mis $FOUT.1stSNP.1stSamp_QC.mis.frqdiff

# Remove A/T C/G SNPs
awk '(($5=="C"&&$6=="G")||($5=="G"&&$6=="C")||($5=="A"&&$6=="T")||($5=="T"&&$6=="A")) {print $2} ' $FOUT.1stSNP.1stSamp_QC.mis.bim > AT_CG_SNPS.txt
plink --bfile $FOUT.1stSNP.1stSamp_QC.mis --exclude AT_CG_SNPS.txt --make-bed --out $FOUT.1stSNP.1stSamp_QC.mis_AT_CG_removed
# 2508 AT_CG_SNPS.txt
bash ${SCRIPTS}/CheckAlleleFreq.sh $FOUT.1stSNP.1stSamp_QC.mis_AT_CG_removed $FOUT.1stSNP.1stSamp_QC.mis.frqdiff_AT_CG_removed

# In R
name = '1K1K'
dat = read.table(paste0(name, ".1stSNP.1stSamp_QC.mis.frqdiff_AT_CG_removed.plot"))
png(paste0(name, ".1stSNP.1stSamp_QC.mis.frqdiff_AT_CG_removed.1KGp3.png"))
plot(dat$V2, dat$V3, pch=20, col="#185fad70", xlab="Freq in Cohort", ylab="Freq in 1KG")
abline(0, 1, col="#e15239", lwd=2, lty=2)
abline(0.3, 1, col="#e15239", lwd=2, lty=2)
abline(-0.3, 1, col="#e15239", lwd=2, lty=2)
dev.off()

# Exclude
rem = dat[abs(dat$V2-dat$V3)>0.3,]
print(nrow(rem))
    
# excluded 2810 SNPs
    
write.table(rem$V1, paste0(name, ".1stSNP.1stSamp_QC.mis.frqdiff_over.0.25.snp.tmp"), quote=F,col.names=F,row.names=F)

# note that the file is labeled "0.25" but should say "0.3"
###########

LANG=en_EN

# Make MAF1% file
awk '{printf("%s_%s\t%s\n",$1,$4,$2)}' $FOUT.1stSNP.1stSamp_QC.mis_AT_CG_removed.bim | sort -k1,1 | join -1 1 -2 1 - $FOUT.1stSNP.1stSamp_QC.mis.frqdiff_over.0.25.snp.tmp | awk '{print $2}' > $FOUT.1stSNP.1stSamp_QC.mis.frqdiff_over.0.25.snp
plink --bfile $FOUT.1stSNP.1stSamp_QC.mis_AT_CG_removed --exclude $FOUT.1stSNP.1stSamp_QC.mis.frqdiff_over.0.25.snp --make-bed --maf 0.01 --hwe 1e-10 --mind 0.01 --geno 0.02 --out $FOUT.final.maf0.01
# --exclude: 490703 variants remaining.
# 0 variants removed due to Hardy-Weinberg exact test.
# 3232 variants removed due to MAF threshold(s)
# 487471 variants and 972 people pass filters and QC.

# Make MAF5% file
awk '{printf("%s_%s\t%s\n",$1,$4,$2)}' $FOUT.1stSNP.1stSamp_QC.mis_AT_CG_removed.bim | sort -k1,1 | join -1 1 -2 1 - $FOUT.1stSNP.1stSamp_QC.mis.frqdiff_over.0.25.snp.tmp | awk '{print $2}' > $FOUT.1stSNP.1stSamp_QC.mis.frqdiff_over.0.25.snp
plink --bfile $FOUT.1stSNP.1stSamp_QC.mis_AT_CG_removed --exclude $FOUT.1stSNP.1stSamp_QC.mis.frqdiff_over.0.25.snp --make-bed --maf 0.05 --hwe 1e-10 --mind 0.01 --geno 0.02 --out $FOUT.final
# --exclude: 490703 variants remaining
# 0 variants removed due to Hardy-Weinberg exact test.
# 212539 variants removed due to MAF threshold(s)
# 278164 variants and 972 people pass filters and QC.

# Ready for HLA imputation: /data/srlab/external-data/1K1K_raw/geno/hg19/qc/1K1K.final.maf0.01.MHC
# Fix sample names in .fam files
cp $FOUT.final.maf0.01.fam $FOUT.final.maf0.01.fam.orig
awk '{ print $1" "$1"_"$2" "$3" "$4" "$5" "$6}' $FOUT.final.maf0.01.fam.orig > $FOUT.final.maf0.01.fam
cp $FOUT.final.fam $FOUT.final.fam.orig
awk '{ print $1" "$1"_"$2" "$3" "$4" "$5" "$6}' $FOUT.final.fam.orig > $FOUT.final.fam