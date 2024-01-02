#!/bin/sh


cd ..



### Step 1: quality control
# fastp installation
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
# run fastp
cd raw_fastq
../fastp -i F1_1_R1.fq.gz -I F1_1_R2.fq.gz -o F1_out_R1.fq.gz -O F1_out_R2.fq.gz
cd ..

### Step 2: alignment with STAR
module add STAR/2.7.9a
# reference build
cd ref
STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles chr2L.fa --sjdbGTFfile chr2L.gtf --runThreadN 1 --genomeSAindexNbases 11
cd ..
# alignment
mkdir bam
cd bam
STAR --genomeDir ../ref --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --alignSJoverhangMin 8 --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c --runThreadN 10 --readFilesIn ../raw_fastq/F1_out_R1.fq.gz ../raw_fastq/F1_out_R2.fq.gz --outFileNamePrefix F1_1
cd ..


### Step 3: gene read count with featureCounts
# featureCounts installation
wget https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-Linux-x86_64.tar.gz
tar -xvzf subread-2.0.6-Linux-x86_64.tar.gz
# run featureCounts
cd bam
../subread-2.0.6-Linux-x86_64/bin/featureCounts -a ../ref/chr2L.gtf -F 'GTF' -t exon -g gene_id -p --countReadPairs -B -Q 20 -o F1_1_count.txt F1_1Aligned.sortedByCoord.out.bam

