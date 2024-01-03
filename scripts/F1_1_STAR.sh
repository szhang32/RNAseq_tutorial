#!/bin/sh
#BSUB -J F1_1_STAR
#BSUB -oo F1_1_STAR.o
#BSUB -eo F1_1_STAR.e
#BSUB -B
#BSUB -N
#BSUB -M 10000
module add STAR/2.7.9a

cd ../bam
STAR --genomeDir ../ref --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --alignSJoverhangMin 8 --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c --runThreadN 10 --readFilesIn ../raw_fastq/F1_out_R1.fq.gz ../raw_fastq/F1_out_R2.fq.gz --outFileNamePrefix F1_1
