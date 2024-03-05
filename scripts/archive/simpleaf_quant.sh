#!/bin/bash
#$ -N simpleafquant
#$ -cwd
#$ -pe omp 15


export ALEVIN_FRY_HOME='/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/alevinfry'
simpleaf set-paths

data_dir="/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/data"
output_dir="/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/alevinfry"
reference_dir="/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/reference"


simpleaf quant \
-c 10xv3 \
-o $output_dir \
-i $output_dir/index \
-1 $data_dir/$1 \
-2 $data_dir/$2  \
-u -r cr-like -m $output_dir/index/t2g_3col.tsv

