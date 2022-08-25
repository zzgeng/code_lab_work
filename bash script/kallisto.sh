#!/bin/bash -l
#SBATCH -p compute
#SBATCH -e kallisto.%j.err
#SBATCH -o kallisto.%j.out

module load kallisto

for x in ZG01 ZG02 ZG03 ZG04 ZG05 ZG06 ZG07 ZG08 ZG09;

do 
kallisto quant -i ~/Kallisto_index/homo_sapiens/transcriptome.idx ${x}/${x}*fastq* -o kallisto_results/${x} 

done
