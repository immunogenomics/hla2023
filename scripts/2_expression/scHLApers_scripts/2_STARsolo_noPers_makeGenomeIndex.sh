#!/bin/bash
## Amber and Joyce
## 19 March 2022
## Run STARsolo (standard pipeline without personalization) genome file generation

module load samtools

# Inputs
genomeDir="/data/srlab2/ashen/hla/scHLA_STARsolo/data/GRCh38_refs/" # store index
numThreads=16

##### Run Non-Personalized Pipeline #####
# Reference genome
GRCh38_genome="/data/srlab2/ashen/hla/scHLA_STARsolo/data/GRCh38_refs/GRCh38.primary_assembly.genome.fa"
GRCh38_annot="/data/srlab1/amber_joyce/filtered_gtf/gencode.v38.annotation.filtered.gtf"

# Generate Genome Index (only has to be run once)
CMD="STAR --runThreadN $numThreads --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $GRCh38_genome --sjdbGTFfile $GRCh38_annot"
$CMD
