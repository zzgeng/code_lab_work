#!/bin/bash -l
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -e callPeaks.%j.err
#SBATCH -o callPeaks.%j.out

mkdir callPeaks	

module load macs2
for file in *.sorted.bam;
do 

macs2 callpeak -t ${file} -f BAM -g mm -n ${file} --outdir callPeaks/

done
