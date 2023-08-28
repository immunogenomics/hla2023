#!/bin/bash
## Amber and Joyce
## 15 August 2022

## Script that runs both the personalized and non-personalized pipeline for a Randolph2021 sample,
## as well as long and short reads separately. Input is the sample name. Script outputs all intermediary files
## and will automatically rerun any job that fails to generate the appropriate output files. 

module load samtools

# Inputs
sample=$1 # sample name
UMIlen=10 # 12 for v3, 10 for v2 (default), Randolph2021 samples are v2

# Set up input files and parameters
readsDir="/data/srlab1/amber_joyce/scHLA/Randolph2021_NewPanel/demultiplexed_bams/"
whitelist="/apps/lib-osver/cellranger/1.3.1/cellranger-cs/1.3.1/lib/python/cellranger/barcodes/737K-august-2016.txt" # v2 whitelist
numThreads=4
soloType="CB_UMI_Simple" # type of CB/UMI used
input_bam=${readsDir}${sample}_shuffled.bam
#samtools collate ${readsDir}${sample}.bam ${readsDir}${sample}_shuffled # already done

# Path to script to make genesXcells
make_exp_script="/data/srlab1/jkang/hla/schla/scripts/2_expression/starsolo_to_genesXcells.R"
cell_meta="/data/srlab1/jkang/hla/data/combined_AMP_Randolph_Smillie/cell_meta_Randolph_completeHLA.csv"

##########################
# Run scHLApers Pipeline
##########################
echo 'Run scHLApers'

label="GeneFull_Ex50pAS"
dir="/data/srlab2/jkang/scHLA/personalized_final/Randolph2021_NewPanel/"
out="${dir}STARsolo_results/${sample}_${label}/${sample}_${label}_" # directory to store results
genomeDir="${dir}${sample}_${label}_index/" # store index

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

    perGenome="${dir}personalized_references/genomes/${sample%%_*}_genome.fa" # Need to remove the NI/flu from sample string
    perAnnot="${dir}personalized_references/nonunique_annotations/${sample%%_*}_annotation.gtf"

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
    CMD="/PHShome/jbk37/tools/STAR-2.7.10a/source/STAR --genomeDir $genomeDir --readFilesIn $input_bam \
          --runThreadN $numThreads --runDirPerm All_RWX --soloUMIlen $UMIlen \
          --outFileNamePrefix $out --soloType $soloType --soloCBwhitelist $whitelist \
          --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR \
          --soloFeatures GeneFull_Ex50pAS \
          --soloBarcodeReadLength 0 \
          --soloMultiMappers EM --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
          --readFilesType SAM SE --readFilesCommand samtools view -F 0x100 \
          --soloInputSAMattrBarcodeSeq CR UR --soloInputSAMattrBarcodeQual CY UY \
          --readFilesSAMattrKeep None \
          --outSAMattributes NH HI nM AS CR UR GX GN CB UB \
          --outTmpDir ${out}_STARtmp" # added temp dir # do not check barcode read length (29 bp)
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

# Delete intermediate files to save space - don't delete personalized genome until after running long/short below
#rm $combPerGenome
#rm $combPerAnnot
#rm -R $genomeDir

CMD="Rscript ${make_exp_script} $sample ${dir}STARsolo_results/ $label $cell_meta"
echo $CMD
$CMD

##########################
# Run scHLApers Pipeline on LONG reads
##########################
echo 'Run scHLApers on long reads'

# Inputs
label="GeneFull_Ex50pAS_longReads"
dir="/data/srlab2/jkang/scHLA/personalized_final/Randolph2021_NewPanel/"
out="${dir}STARsolo_results/${sample}_${label}/${sample}_${label}_" # directory to store results

readsDir="/data/srlab2/jkang/scHLA/Randolph_long_short_bams/"
input_bam=${readsDir}${sample}_shuffled_long.bam

finished=false
while ! $finished; do
    ##### Run Personalized Pipeline #####
    rm -rf ${dir}STARsolo_results/${sample}_${label} # remove the output directory if it exists
    mkdir -p ${dir}STARsolo_results/${sample}_${label}
    
    # Genome already generated previously

    # Align Reads
    CMD="/PHShome/jbk37/tools/STAR-2.7.10a/source/STAR --genomeDir $genomeDir --readFilesIn $input_bam \
          --runThreadN $numThreads --runDirPerm All_RWX --soloUMIlen $UMIlen \
          --outFileNamePrefix $out --soloType $soloType --soloCBwhitelist $whitelist \
          --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR \
          --soloFeatures GeneFull_Ex50pAS \
          --soloBarcodeReadLength 0 \
          --soloMultiMappers EM --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
          --readFilesType SAM SE --readFilesCommand samtools view -F 0x100 \
          --soloInputSAMattrBarcodeSeq CR UR --soloInputSAMattrBarcodeQual CY UY \
          --readFilesSAMattrKeep None \
          --outSAMattributes NH HI nM AS CR UR GX GN CB UB \
          --outTmpDir ${out}_STARtmp" # added temp dir # do not check barcode read length (29 bp)
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

CMD="Rscript ${make_exp_script} $sample ${dir}STARsolo_results/ $label $cell_meta"
echo $CMD
$CMD

###############
## Run scHLApers Pipeline on SHORT reads
###############
echo 'Run scHLApers on short reads'

# Inputs
label="GeneFull_Ex50pAS_shortReads"
dir="/data/srlab2/jkang/scHLA/personalized_final/Randolph2021_NewPanel/"
out="${dir}STARsolo_results/${sample}_${label}/${sample}_${label}_" # directory to store results

readsDir="/data/srlab2/jkang/scHLA/Randolph_long_short_bams/"
input_bam=${readsDir}${sample}_shuffled_short.bam

finished=false
while ! $finished; do

    ##### Run Personalized Pipeline #####
    rm -rf ${dir}STARsolo_results/${sample}_${label} # remove the output directory if it exists
    mkdir -p ${dir}STARsolo_results/${sample}_${label}

    # Align Reads
    CMD="/PHShome/jbk37/tools/STAR-2.7.10a/source/STAR --genomeDir $genomeDir --readFilesIn $input_bam \
          --runThreadN $numThreads --runDirPerm All_RWX --soloUMIlen $UMIlen \
          --outFileNamePrefix $out --soloType $soloType --soloCBwhitelist $whitelist \
          --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR \
          --soloFeatures GeneFull_Ex50pAS \
          --soloBarcodeReadLength 0 \
          --soloMultiMappers EM --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
          --readFilesType SAM SE --readFilesCommand samtools view -F 0x100 \
          --soloInputSAMattrBarcodeSeq CR UR --soloInputSAMattrBarcodeQual CY UY \
          --readFilesSAMattrKeep None \
          --outSAMattributes NH HI nM AS CR UR GX GN CB UB \
          --outTmpDir ${out}_STARtmp" # added temp dir # do not check barcode read length (29 bp)
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


#################################################################################################
###############
## Run standard pipeline
###############
echo 'Run standard'

label="GeneFull_Ex50pAS"
dir="/data/srlab2/jkang/scHLA/noPers_final/Randolph2021_NewPanel/"
genomeDir="/data/srlab2/ashen/hla/scHLA_STARsolo/data/GRCh38_refs/" # store index
out="${dir}STARsolo_results/${sample}_${label}/${sample}_${label}_" # directory to store results

readsDir="/data/srlab1/amber_joyce/scHLA/Randolph2021_NewPanel/demultiplexed_bams/"
input_bam=${readsDir}${sample}_shuffled.bam

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
    CMD="/PHShome/jbk37/tools/STAR-2.7.10a/source/STAR --genomeDir $genomeDir --readFilesIn $input_bam \
          --runThreadN $numThreads --runDirPerm All_RWX --soloUMIlen $UMIlen \
          --outFileNamePrefix $out --soloType $soloType --soloCBwhitelist $whitelist \
          --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR \
          --soloFeatures GeneFull_Ex50pAS \
          --soloBarcodeReadLength 0 \
          --soloMultiMappers EM --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
          --readFilesType SAM SE --readFilesCommand samtools view -F 0x100 \
          --soloInputSAMattrBarcodeSeq CR UR --soloInputSAMattrBarcodeQual CY UY \
          --readFilesSAMattrKeep None \
          --outSAMattributes NH HI nM AS CR UR GX GN CB UB" # do not check barcode read length (29 bp)
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

###############
## Run standard on LONG reads
###############
echo 'Run standard on long reads'

label="GeneFull_Ex50pAS_longReads"
out="${dir}STARsolo_results/${sample}_${label}/${sample}_${label}_" # directory to store results

readsDir="/data/srlab2/jkang/scHLA/Randolph_long_short_bams/"
input_bam=${readsDir}${sample}_shuffled_long.bam

finished=false
while ! $finished; do
    rm -rf ${dir}STARsolo_results/${sample}_${label} # remove the output directory if it exists
    mkdir -p ${dir}STARsolo_results/${sample}_${label}

    # Align Reads
    CMD="/PHShome/jbk37/tools/STAR-2.7.10a/source/STAR --genomeDir $genomeDir --readFilesIn $input_bam \
          --runThreadN $numThreads --runDirPerm All_RWX --soloUMIlen $UMIlen \
          --outFileNamePrefix $out --soloType $soloType --soloCBwhitelist $whitelist \
          --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR \
          --soloFeatures GeneFull_Ex50pAS \
          --soloBarcodeReadLength 0 \
          --soloMultiMappers EM --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
          --readFilesType SAM SE --readFilesCommand samtools view -F 0x100 \
          --soloInputSAMattrBarcodeSeq CR UR --soloInputSAMattrBarcodeQual CY UY \
          --readFilesSAMattrKeep None \
          --outSAMattributes NH HI nM AS CR UR GX GN CB UB" # do not check barcode read length (29 bp)
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


###############
## Run standard on SHORT reads
###############
echo 'Run standard on short reads'

label="GeneFull_Ex50pAS_shortReads"
out="${dir}STARsolo_results/${sample}_${label}/${sample}_${label}_" # directory to store results

readsDir="/data/srlab2/jkang/scHLA/Randolph_long_short_bams/"
input_bam=${readsDir}${sample}_shuffled_short.bam

finished=false
while ! $finished; do
    rm -rf ${dir}STARsolo_results/${sample}_${label} # remove the output directory if it exists
    mkdir -p ${dir}STARsolo_results/${sample}_${label}

    # Align Reads
    CMD="/PHShome/jbk37/tools/STAR-2.7.10a/source/STAR --genomeDir $genomeDir --readFilesIn $input_bam \
          --runThreadN $numThreads --runDirPerm All_RWX --soloUMIlen $UMIlen \
          --outFileNamePrefix $out --soloType $soloType --soloCBwhitelist $whitelist \
          --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR \
          --soloFeatures GeneFull_Ex50pAS \
          --soloBarcodeReadLength 0 \
          --soloMultiMappers EM --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
          --readFilesType SAM SE --readFilesCommand samtools view -F 0x100 \
          --soloInputSAMattrBarcodeSeq CR UR --soloInputSAMattrBarcodeQual CY UY \
          --readFilesSAMattrKeep None \
          --outSAMattributes NH HI nM AS CR UR GX GN CB UB" # do not check barcode read length (29 bp)
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