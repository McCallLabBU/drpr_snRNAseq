#!/bin/bash
#$ -N cellrangercount
#$ -cwd
#$ -pe omp 36

module load bcl2fastq/2.20                                                                                                                                   
module load cellranger/7.2.0

cd /projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis

cellranger count --id=$1 \
   --fastqs=/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/data \
   --sample=$2 \
   --transcriptome=/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/reference/BDGP6.32.109 \
   --include-introns true \
   --expect-cells 10000