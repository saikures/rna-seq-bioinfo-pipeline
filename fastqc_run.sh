#!/bin/bash
#PBS -N fastqc_post_maj_strin
#PBS -l select=1:ncpus=126:mem=100gb
#PBS -l walltime=24:00:00
#PBS -o fastqc_post_maj_strin.log
#PBS -j oe

export OMP_NUM_THREADS=126

OUTDIR=/storage/brno12-cerit/home/saikures/saikures_main/RNASeq_major/maj_dataset/trim_maj_stringent/fastqc_trim_strin_maj

cd "$PBS_O_WORKDIR"
mkdir "$OUTDIR"

source /storage/brno12-cerit/home/saikures/miniconda3/bin/activate
source ~/.bashrc

conda activate fastqc_env

fastqc M-*.fastq.gz -t 126 -o "$OUTDIR"

