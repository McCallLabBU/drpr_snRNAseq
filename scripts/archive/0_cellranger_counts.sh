#!/bin/bash
#$ -N cellrangercount
#$ -cwd
#$ -pe omp 36

module load bcl2fastq/2.20                                                                                                                                   
module load cellranger/6.1.2

cd /projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis

cellranger count --id=drpr_w1118_42d_SNrnaseq \
   --fastqs=/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/data \
   --sample=Drprnull42D_CKDL230005442-1A_HWKMGDSX5,W1118_42D_CKDL230005441-1A_HWKMGDSX5 \
   --transcriptome=/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/reference/BDGP6.32.109