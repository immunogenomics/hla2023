#!/bin/bash
## Amber and Joyce
## 20 May 2022
## Run STARsolo (without personalization) for Smillie cohort

module load samtools

# Inputs
sample=$1 # sample name
UMIlen=$2 # Smillie samples are all v1 or v2 (10 bp)
version=$3 # either "v1" or "v2"

label="GeneFull_Ex50pAS"
dir="/data/srlab2/jkang/scHLA/noPers_final/Smillie2019_NewPanel/"
genomeDir="/data/srlab2/ashen/hla/scHLA_STARsolo/data/GRCh38_refs/" # store index
out="${dir}STARsolo_results/${sample}_${label}/${sample}_${label}_" # directory to store results
mkdir -p ${dir}STARsolo_results/${sample}_${label}

readsDir="/data/srlab1/jkang/hla/data/Smillie2019/gcloud/bams/"

if [[ $version == "v1" ]]; then
   echo "using v1 format"
   whitelist="/apps/lib-osver/cellranger/1.3.1/cellranger-cs/1.3.1/lib/python/cellranger/barcodes/737K-april-2014_rc.txt" # v1 3' whitelist
   CBlen=14
elif [[ $version == "v2" ]]; then
   echo "using v2 format"
   whitelist="/apps/lib-osver/cellranger/1.3.1/cellranger-cs/1.3.1/lib/python/cellranger/barcodes/737K-august-2016.txt" # v2 3' whitelist
   CBlen=16
else
   echo "invalid 10x version $version"
fi

numThreads=4
soloType="CB_UMI_Simple" # type of CB/UMI used
genomeDir="/data/srlab2/ashen/hla/scHLA_STARsolo/data/GRCh38_refs/" # store index

##### Run Non-Personalized Pipeline #####
for cond in EpiA EpiB LPA LPB
do
   cond_out="${out}${cond}."
   input_bam=${readsDir}${sample}.${cond}_shuffled.bam
   samtools collate -@ 4 ${readsDir}${sample}.${cond}.bam ${readsDir}${sample}.${cond}_shuffled
   
   # Reference genome
   #GRCh38_genome="/data/srlab2/ashen/hla/scHLA_STARsolo/data/GRCh38_refs/GRCh38.primary_assembly.genome.fa"
   #GRCh38_annot="/data/srlab1/amber_joyce/filtered_gtf/gencode.v38.annotation.filtered.gtf"
   
   # Generate Genome Index (only has to be run once, comment out if already run)
   #CMD="STAR --runThreadN $numThreads --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $GRCh38_genome --sjdbGTFfile $GRCh38_annot"
   #$CMD
   
   # Align Reads
    CMD="/PHShome/jbk37/tools/STAR-2.7.10a/source/STAR --genomeDir $genomeDir --readFilesIn $input_bam \
          --runThreadN $numThreads --runDirPerm All_RWX --soloUMIlen $UMIlen --soloCBlen $CBlen \
          --outFileNamePrefix $cond_out --soloType $soloType --soloCBwhitelist $whitelist \
          --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR \
          --soloFeatures GeneFull_Ex50pAS \
          --soloMultiMappers EM --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
          --readFilesType SAM SE --readFilesCommand samtools view -F 0x100 \
          --soloInputSAMattrBarcodeSeq CR UR --soloInputSAMattrBarcodeQual CY UY \
          --readFilesSAMattrKeep None \
          --outSAMattributes NH HI nM AS CR UR GX GN CB UB"
    echo $CMD
    $CMD

    # Sort the output bam file and save HLA region and unmapped reads
    bam="${cond_out}Aligned.sortedByCoord.out.bam"
    CMD="samtools index -@ 6 ${bam}"
    echo $CMD
    $CMD

    bam_HLA="${cond_out}HLA.bam"
    bam_unmapped="${cond_out}unmapped.bam"
   
    CMD="samtools view -b ${bam} chr6:28000000-34000000 -o ${bam_HLA}"
    $CMD
    CMD="samtools view -b -f 4 ${bam} -o ${bam_unmapped}"
    $CMD
    CMD="samtools index -@ 6 ${bam_HLA}"
    $CMD

    # Remove original output bam file to save space
    rm $bam
    rm "${bam}.bai"
done
