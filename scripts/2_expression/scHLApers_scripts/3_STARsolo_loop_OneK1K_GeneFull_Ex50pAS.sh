#!/bin/bash
## Amber and Joyce
## 15 August 2022

## Script that runs both the personalized and non-personalized pipeline for a OneK1K sample.
## Input is the sample name and path to STAR executable. Script outputs all intermediary files
## and will automatically rerun any job that fails to generate the appropriate output files.

module load samtools

# Inputs
sample=$1 # sample name
STAR_executable=$2 # e.g. /PHShome/jbk37/tools/STAR-2.7.10a/source/STAR
UMIlen=10 # 10 for v2
label="GeneFull_Ex50pAS" # Label for the output directory

# Set up input files and parameters
readsDir="/data/srlab/external-data/1K1K_raw/demultiplexed_bams/"
whitelist="/apps/lib-osver/cellranger/1.3.1/cellranger-cs/1.3.1/lib/python/cellranger/barcodes/737K-august-2016.txt" # v2 whitelist
numThreads=4
soloType="CB_UMI_Simple" # type of CB/UMI used
input_bam=${readsDir}${sample}_shuffled.bam
samtools collate ${readsDir}${sample}.bam ${readsDir}${sample}_shuffled # Shuffle the input file

# Path to script to make genesXcells
make_exp_script="/data/srlab1/jkang/hla/schla/scripts/2_expression/starsolo_to_genesXcells.R"
cell_meta="/data/srlab1/jkang/hla/data/combined_AMP_Randolph_Smillie/cell_meta_OneK1K_completeHLA.csv"

###########################
### RUN scHLApers PIPELINE
###########################

# Set up output directory
dir="/data/srlab2/jkang/scHLA/personalized_final/OneK1K_NewPanel/"
out="${dir}STARsolo_results/${sample}_${label}/${sample}_${label}_" # directory to store results
genomeDir="${dir}${sample}_${label}_index/" # where to store genome index

finished=false
while ! $finished; do
    ##### Run Personalized Pipeline #####
    rm -rf ${dir}STARsolo_results/${sample}_${label} # remove the output directory if it exists
    mkdir -p ${dir}STARsolo_results/${sample}_${label}
    
    rm -rf $genomeDir # remove the genome directory if it exists
    mkdir $genomeDir
    
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
    CMD="$STAR_executable --runThreadN $numThreads \
        --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $combPerGenome \
        --sjdbGTFfile $combPerAnnot --outTmpDir ${out}_STARtmp_genome_generate"
    echo $CMD
    $CMD

    # Align Reads
    CMD="$STAR_executable --genomeDir $genomeDir --readFilesIn $input_bam \
          --runThreadN $numThreads --runDirPerm All_RWX --soloUMIlen $UMIlen \
          --outFileNamePrefix $out --soloType $soloType --soloCBwhitelist $whitelist \
          --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR \
          --soloFeatures GeneFull_Ex50pAS \
          --soloMultiMappers EM --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
          --readFilesType SAM SE --readFilesCommand samtools view -F 0x100 \
          --soloInputSAMattrBarcodeSeq CR UR --soloInputSAMattrBarcodeQual CY UY \
          --readFilesSAMattrKeep None \
          --outSAMattributes NH HI nM AS CR UR GX GN CB UB --outTmpDir ${out}_STARtmp" # added temp dir
    echo $CMD
    $CMD
    
    ## Check whether the run generated output
    if [ -f "${out}Solo.out/GeneFull_Ex50pAS/raw/UniqueAndMult-EM.mtx" ] && [ -f "${out}Aligned.sortedByCoord.out.bam" ]; then
        echo "File exists, setting finished=true"
        finished=true
    else
        echo "Error: File does not exist, rerunning STARsolo personalized pipeline"
    fi
done

# Remove shuffled version of input bam (save space)
#rm $input_bam # Keep shuffled input bam for standard pipeline below

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

###################################################################################################################
###########################
### RUN Standard PIPELINE
###########################

dir="/data/srlab2/jkang/scHLA/noPers_final/OneK1K_NewPanel/"
genomeDir="/data/srlab2/ashen/hla/scHLA_STARsolo/data/GRCh38_refs/" # store index
out="${dir}STARsolo_results/${sample}_${label}/${sample}_${label}_" # directory to store results

# Input bam same as above (already shuffled)

##### Run Non-Personalized Pipeline #####
# Reference genome
#GRCh38_genome="/data/srlab2/ashen/hla/scHLA_STARsolo/data/GRCh38_refs/GRCh38.primary_assembly.genome.fa"
#GRCh38_annot="/data/srlab1/amber_joyce/filtered_gtf/gencode.v38.annotation.filtered.gtf"

# Generate Genome Index (only has to be run once, comment out if already run)
#CMD="STAR --runThreadN $numThreads --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $GRCh38_genome --sjdbGTFfile $GRCh38_annot"
#$CMD

finished=false
while ! $finished; do
    rm -rf ${dir}STARsolo_results/${sample}_${label} # remove the output directory if it exists
    mkdir -p ${dir}STARsolo_results/${sample}_${label}
    
    # Align Reads
    CMD="$STAR_executable --genomeDir $genomeDir --readFilesIn $input_bam \
          --runThreadN $numThreads --runDirPerm All_RWX --soloUMIlen $UMIlen \
          --outFileNamePrefix $out --soloType $soloType --soloCBwhitelist $whitelist \
          --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR \
          --soloFeatures GeneFull_Ex50pAS \
          --soloMultiMappers EM --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
          --readFilesType SAM SE --readFilesCommand samtools view -F 0x100 \
          --soloInputSAMattrBarcodeSeq CR UR --soloInputSAMattrBarcodeQual CY UY \
          --readFilesSAMattrKeep None \
          --outSAMattributes NH HI nM AS CR UR GX GN CB UB"
    echo $CMD
    $CMD
    
    ## Check whether the run generated output
    if [ -f "${out}Solo.out/GeneFull_Ex50pAS/raw/UniqueAndMult-EM.mtx" ] && [ -f "${out}Aligned.sortedByCoord.out.bam" ]; then
        echo "File exists, setting finished=true"
        finished=true
    else
        echo "Error: File does not exist, rerunning STARsolo personalized pipeline"
    fi
done

# Remove shuffled version of input bam (save space)
rm $input_bam

# Sort the output bam file and save HLA region and unmapped reads
bam="${out}Aligned.sortedByCoord.out.bam"
CMD="samtools index -@ 6 ${bam}"
echo $CMD
$CMD

bam_HLA="${out}HLA.bam"
bam_unmapped="${out}unmapped.bam"
   
CMD="samtools view -b ${bam} chr6:28000000-34000000 -o ${bam_HLA}"
$CMD
CMD="samtools view -b -f 4 ${bam} -o ${bam_unmapped}"
$CMD
CMD="samtools index -@ 6 ${bam_HLA}"
$CMD

# Remove original output bam file to save space
rm $bam
rm "${bam}.bai"

CMD="Rscript ${make_exp_script} $sample ${dir}STARsolo_results/ $label $cell_meta"
echo $CMD
$CMD

echo "Finished successfully."