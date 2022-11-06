#!/bin/bash

#define the paths of the raw sequence data and the T. congolense genome sequence as variables 
FASTQ_PATH=/localdisk/data/BPSM/AY21/fastq
GENOME_PATH=/localdisk/data/BPSM/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz

#create new directories to store counts data outputs, intermediate alignment .bam files, and the genome_index
mkdir counts_data
mkdir bam_files
mkdir genome_index

#creates a new genome index for the Trypanasoma genome using bowtie-build
echo "Creating a genome index..."
bowtie2-build $GENOME_PATH genome_index/Tcongo_genome_index --threads 16 &&
echo "Genome index built!"

#loop through each fastq file to pick out the paried end reads and use bowtie2 to align them to the T. congolense genome and pipe the .sam output to samtools view converting to compressed binary .bam format
for i in ${FASTQ_PATH}/*1.fq.gz 
do
    SAMPLE=$(basename $i _1.fq.gz); 
    echo "Aligning and counting ${SAMPLE}"
    bowtie2 -p 16 -x genome_index/Tcongo_genome_index -1 ${FASTQ_PATH}/${SAMPLE}_1.fq.gz -2 ${FASTQ_PATH}/${SAMPLE}_2.fq.gz | samtools view -bS > bam_files/${SAMPLE}.bam &&

#still within loop, samtools sort sorts the .bam files to produce .sorted.bam files 
    samtools sort bam_files/${SAMPLE}.bam -o bam_files/${SAMPLE}_sorted.bam &&
 
#samtools index takes sorted.bam files and creates index files for fast access
    samtools index bam_files/${SAMPLE}_sorted.bam && 

#the bedtools coverage command takes the sorted.bam files and counts the number of alignments corresponding to gene regions in the genome using the reference .bed file that contains infomration on gene locations
#the resulting counts files are stored in the counts_data directory
    bedtools coverage -counts -a /localdisk/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed -b bam_files/${SAMPLE}_sorted.bam -F 0.1 > counts_data/${SAMPLE}_counts.txt
done &&

echo "Counts data generated!"
