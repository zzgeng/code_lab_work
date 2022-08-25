#!/bin/bash -l
#SBATCH -N1
#SBATCH --partition=compute
#SBATCH -o chip_TAM_H2A.out
#SBATCH -e chip_TAM_H2A.err
#SBATCH --mail-user=zgeng@pennstatehealth.psu.edu
#SBATCH --mail-type=ALL

module load sratoolkit
module load bwa
module load samtools
module load deeptools

fastq-dump SRX4650800
mv SRX4650800.fastq P2_TAM_H2Aub_1.fastq

fastq-dump SRX4650801
mv SRX4650801.fastq P24_TAM_H2Aub_2.fastq

fastq-dump SRX4650802
mv SRX4650802.fastq P24_TAM_H2Aub_3.fastq


bwa mem -t 16 ~/bwa_index/mm10/genome.fa P2_TAM_H2Aub_1.fastq | samtools view -bS - > P2TAM_H2A_1.bam
bwa mem -t 16 ~/bwa_index/mm10/genome.fa P2_TAM_H2Aub_2.fastq | samtools view -bS - > P2TAM_H2A_2.bam
bwa mem -t 16 ~/bwa_index/mm10/genome.fa P2_TAM_H2Aub_3.fastq | samtools view -bS - > P2TAM_H2A_3.bam

samtools merge P2TAM_H2A.bam P2TAM_H2A_1.bam P2TAM_H2A_2.bam P2TAM_H2A_3.bam
samtools sort P2TAM_H2A.bam > P2TAM_H2A.sorted.bam 
samtools index P2TAM_H2A.sorted.bam

bamCoverage -b P2TAM_H2A.sorted.bam -o P2TAM_H2A.BPM.bw --normalizeUsing BPM --effectiveGenomeSize 2308125349
