#!/bin/bash
#$ -N cellranger
#$ -cwd
#$ -pe omp 36

module load bcl2fastq/2.20                                                                                                                                   
module load cellranger/6.1.2

cd /projectnb/mccall/sbandyadka/drpr42d_snrnaseq/reference/

#Remove previous runs
#rm -rf /projectnb/mccall/sbandyadka/drpr42d_snrnaseq/reference/BDGP6.32.109


cellranger mkgtf /projectnb/mccall/sbandyadka/drpr42d_snrnaseq/reference/Drosophila_melanogaster.BDGP6.32.109.gtf /projectnb/mccall/sbandyadka/drpr42d_snrnaseq/reference/Drosophila_melanogaster.BDGP6.32.109.filtered.gtf \
   --attribute=gene_biotype:protein_coding \
   --attribute=gene_biotype:ncRNA 

cellranger mkref --genome=BDGP6.32.109 \
    --fasta=/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/reference/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa \
    --genes=/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/reference/Drosophila_melanogaster.BDGP6.32.109.filtered.gtf \
    --nthreads=36