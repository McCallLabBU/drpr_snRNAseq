#!/bin/bash
#$ -N simpleafindex
#$ -cwd
#$ -pe omp 15

#module load miniconda
#conda activate af

export ALEVIN_FRY_HOME='/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/alevinfry'


simpleaf set-paths


output_dir="/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/alevinfry"
reference_dir="/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/reference"

simpleaf index \
-o $output_dir \
-f $reference_dir/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa \
-g $reference_dir/Drosophila_melanogaster.BDGP6.32.109.gtf \
-r 150 \
-t 10 \
--use-piscem