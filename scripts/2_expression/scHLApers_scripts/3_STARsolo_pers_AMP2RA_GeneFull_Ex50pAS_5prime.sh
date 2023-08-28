#!/bin/bash
## Amber and Joyce
## 5 July 2022

module load samtools

# Inputs
sample=$1 # sample name
UMIlen=$2 # 12 for v3 format, 10 for v2 format (use v3 whitelist for all samples)
label="GeneFull_Ex50pAS_5prime"
dir="/data/srlab2/jkang/scHLA/personalized_final/AMP2RA_NewPanel/"
genomeDir="${dir}${sample}_${label}_index/" # store index
rm -rf $genomeDir # remove the genome directory if it exists
mkdir $genomeDir

out="${dir}STARsolo_results/${sample}_${label}/${sample}_${label}_" # directory to store results
rm -rf ${dir}STARsolo_results/${sample}_${label} # remove the output directory if it exists
mkdir ${dir}STARsolo_results/${sample}_${label}

readsDir="/data/srlab2/jkang/scHLA/AMP_RA_repertoire/gex/fastq/inputfq/"
whitelist="/apps/lib-osver/cellranger/1.3.1/cellranger-cs/1.3.1/lib/python/cellranger/barcodes/737K-august-2016.txt" # v2 whitelist
numThreads=4
soloType="CB_UMI_Simple" # type of CB/UMI used
fastq2=${readsDir}${sample}_R2.fastq.gz
fastq1=${readsDir}${sample}_R1.fastq.gz

# Make genesXcells paths
make_exp_script="/data/srlab1/jkang/hla/schla/scripts/2_expression/starsolo_to_genesXcells.R"
cell_meta="/data/srlab1/jkang/hla/longShort_5prime/AMP5prime_metadata_Syn.csv"

##### Run Personalized Pipeline #####
# Combine personalized genome with reference (GRCh38)
GRCh38_genome="/data/srlab2/ashen/hla/scHLA_STARsolo/data/GRCh38_refs/GRCh38.primary_assembly.genome_maskedHLAClassical.fa"
GRCh38_annot="/data/srlab1/amber_joyce/filtered_gtf/gencode.v38.annotation.filtered.masked.gtf"

perGenome="${dir}personalized_references/genomes/${sample}_genome.fa"
perAnnot="${dir}personalized_references/nonunique_annotations/${sample}_annotation.gtf"

combPerGenome="${dir}STARsolo_results/${sample}_${label}_genome.fa"
combPerAnnot="${dir}STARsolo_results/${sample}_${label}_annot.gtf"

cat $GRCh38_genome $perGenome > $combPerGenome
cat $GRCh38_annot $perAnnot > $combPerAnnot

# Generate Genome Index
CMD="/PHShome/jbk37/tools/STAR-2.7.10a/source/STAR --runThreadN $numThreads \
        --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $combPerGenome \
        --sjdbGTFfile $combPerAnnot --outTmpDir ${out}_STARtmp_genome_generate"
$CMD

# Align Reads
CMD="/PHShome/jbk37/tools/STAR-2.7.10a/source/STAR --genomeDir $genomeDir --readFilesIn $fastq2 $fastq1 --runThreadN $numThreads --runDirPerm All_RWX --soloUMIlen $UMIlen \
          --outFileNamePrefix $out --soloType $soloType --soloCBwhitelist $whitelist \
          --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR \
          --soloFeatures GeneFull_Ex50pAS \
          --soloMultiMappers EM --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
          --readFilesCommand zcat \
          --outSAMattributes NH HI nM AS CR UR GX GN CB UB --outTmpDir ${out}_STARtmp" # added temp dir
echo $CMD
$CMD

# Sort the output bam file and save HLA region and unmapped reads
bam="${out}Aligned.sortedByCoord.out.bam"
CMD="samtools index -@ 6 ${bam}"
echo $CMD
$CMD

pers_chrs=$(samtools idxstats ${bam} | cut -f 1 | grep "IMGT_") # extract names from bam file
pers_chrs=$(echo $pers_chrs | sed -e 's/ /: /g') # add ":" after every space
pers_chrs="${pers_chrs}:" # add ":" after final chromosome name

bam_HLA="${out}HLA.bam"
bam_unmapped="${out}unmapped.bam"
   
CMD="samtools view -b ${bam} chr6:28000000-34000000 ${pers_chrs} -o ${bam_HLA}"
$CMD
CMD="samtools view -b -f 4 ${bam} -o ${bam_unmapped}"
$CMD
CMD="samtools index -@ 6 ${bam_HLA}"
$CMD

# Remove original output bam file to save space
rm $bam
rm "${bam}.bai"

# Delete intermediate files to save space
rm $combPerGenome
rm $combPerAnnot
rm -R $genomeDir

CMD="Rscript ${make_exp_script} $sample ${dir}STARsolo_results/ $label $cell_meta"
echo $CMD
$CMD
