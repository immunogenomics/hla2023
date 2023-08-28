#!/bin/bash
## Amber and Joyce
## 17 May 2022
## Run STARsolo (with personalization) for Smillie cohort

# Inputs
sample=$1
UMIlen=$2 # 12 for v3, 10 for v2, 10 for v2; for this cohort should always be 10
version=$3 # either "v1" or "v2" 
label="GeneFull_Ex50pAS"

module load samtools

dir="/data/srlab2/jkang/scHLA/personalized_final/Smillie2019_NewPanel/"
genomeDir="${dir}${sample}_${label}_index/" # store index
rm -rf $genomeDir # remove the genome directory if it exists
mkdir $genomeDir

out="${dir}STARsolo_results/${sample}_${label}/${sample}_${label}_" # directory to store results
rm -rf ${dir}STARsolo_results/${sample}_${label} # remove the output directory if it exists
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

##### Run Personalized Pipeline #####
# Combine personalized genome with reference (GRCh38)
GRCh38_genome="/data/srlab2/ashen/hla/scHLA_STARsolo/data/GRCh38_refs/GRCh38.primary_assembly.genome_maskedHLAClassical.fa"
GRCh38_annot="/data/srlab1/amber_joyce/filtered_gtf/gencode.v38.annotation.filtered.masked.gtf"

perGenome="${dir}personalized_references/genomes/${sample}_genome.fa"
perAnnot="${dir}personalized_references/nonunique_annotations/${sample}_annotation.gtf"

combPerGenome="${dir}STARsolo_results/${sample}_nonunique_genome.fa"
combPerAnnot="${dir}STARsolo_results/${sample}_nonunique_annot.gtf"

cat $GRCh38_genome $perGenome > $combPerGenome
cat $GRCh38_annot $perAnnot > $combPerAnnot

# Generate Genome Index
CMD="/PHShome/jbk37/tools/STAR-2.7.10a/source/STAR --runThreadN $numThreads \
        --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $combPerGenome \
        --sjdbGTFfile $combPerAnnot --outTmpDir ${out}_STARtmp_genome_generate"
$CMD

# Align Reads
for cond in EpiA EpiB LPA LPB
do
   cond_out="${out}${cond}."
   input_bam="${readsDir}${sample}.${cond}_shuffled.bam"
   
   CMD="/PHShome/jbk37/tools/STAR-2.7.10a/source/STAR --genomeDir $genomeDir --readFilesIn $input_bam \
          --runThreadN $numThreads --runDirPerm All_RWX --soloUMIlen $UMIlen --soloCBlen $CBlen \
          --outFileNamePrefix $cond_out --soloType $soloType --soloCBwhitelist $whitelist \
          --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR \
          --soloFeatures GeneFull_Ex50pAS \
          --soloMultiMappers EM --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
          --readFilesType SAM SE --readFilesCommand samtools view -F 0x100 \
          --soloInputSAMattrBarcodeSeq CR UR --soloInputSAMattrBarcodeQual CY UY \
          --readFilesSAMattrKeep None \
          --outSAMattributes NH HI nM AS CR UR GX GN CB UB --outTmpDir ${cond_out}_STARtmp" # added temp dir
   echo $CMD
   $CMD

   # Sort the output bam file and save HLA region and unmapped reads
    bam="${cond_out}Aligned.sortedByCoord.out.bam"
    CMD="samtools index -@ 6 ${bam}"
    echo $CMD
    $CMD

    pers_chrs=$(samtools idxstats ${bam} | cut -f 1 | grep "IMGT_") # extract names from bam file
    pers_chrs=$(echo $pers_chrs | sed -e 's/ /: /g') # add ":" after every space
    pers_chrs="${pers_chrs}:" # add ":" after final chromosome name

    bam_HLA="${cond_out}HLA.bam"
    bam_unmapped="${cond_out}unmapped.bam"
   
    CMD="samtools view -b ${bam} chr6:28000000-34000000 ${pers_chrs} -o ${bam_HLA}"
    $CMD
    CMD="samtools view -b -f 4 ${bam} -o ${bam_unmapped}"
    $CMD
    CMD="samtools index -@ 6 ${bam_HLA}"
    $CMD

    # Remove original output bam file to save space
    rm $bam
    rm "${bam}.bai"
done

# Delete intermediate files to save space
rm $combPerGenome
rm $combPerAnnot
rm -R $genomeDir
