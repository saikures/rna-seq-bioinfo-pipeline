#!/bin/bash

# ==============================================================================
# PBS DIRECTIVES
# ==============================================================================
#PBS -N RNA-Seq-Pipeline      # Job name
#PBS -o snakemake_output.log  # Standard output file
#PBS -e snakemake_error.log   # Standard error file
#PBS -l select=1:ncpus=1      # Request one node with one core for Snakemake's master process
#PBS -l walltime=48:00:00     # Time limit hrs:min:sec
#PBS -l mem=4gb               # Memory for the main Snakemake process
#PBS -q long                  # Specify the queue/partition

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================
# This section remains the same as it sets up the environment regardless of the scheduler.
echo "Loading environment modules..."
module load snakemake/latest
module load fastqc/latest
module load trimmomatic/latest
module load trinity/latest
module load bowtie2/latest
module load samtools/latest
module load kallisto/latest
module load r/latest

# ==============================================================================
# RUN SNAKEMAKE
# ==============================================================================
# The core Snakemake command.
# The `snakemake` command itself is submitted to the cluster, which will
# then submit individual jobs for each rule based on the resources requested.
echo "Starting Snakemake pipeline..."
snakemake --cores 8 \
          --latency-wait 60 \
          --config use_scratchdir=true \
          --resources mem_mb=250000 \
          --jobname '{rule}.{jobid}' \
          --use-conda \
          -p

echo "Snakemake pipeline finished."

