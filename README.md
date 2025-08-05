# RNA-Seq Bioinformatics Pipeline

This repository contains scripts and configurations for RNA-Seq analysis including de novo assembly, scaffolding, and visualization.

RNA-Seq De Novo Transcriptomics Pipeline

This document describes a comprehensive and reproducible bioinformatics pipeline for RNA-Seq analysis using de novo transcriptome assembly. The pipeline is orchestrated using Snakemake and is designed to run efficiently on both local machines and High-Performance Computing (HPC) clusters with a PBS scheduler.

The pipeline takes raw FASTQ files, performs quality control and trimming, assembles a transcriptome using Trinity, quantifies gene expression, and generates various downstream analysis plots using R.
1. Pipeline Overview

The pipeline is structured as a series of interdependent rules, with each step building on the previous one.

    Quality Control (QC):

        Initial QC on raw reads using FastQC.

        Aggregation of QC reports with MultiQC.

    Read Trimming:

        Removal of adapter sequences and low-quality bases using Trimmomatic.

    De Novo Assembly:

        Assembly of the trimmed reads into a transcriptome using Trinity.

    Quantification:

        Indexing of the Trinity assembly and quantification of reads using Kallisto.

    Statistical Analysis & Visualization:

        A dedicated R script (rna_analysis.R) performs differential expression analysis with DESeq2.

        Generates various plots, including PCA, MA, Volcano, p-value distribution, gene regulation plots, and heatmaps.

    Visualization of Reads:

        Aligns reads back to the Trinity assembly using Bowtie2 to create BAM files.

        Indexes the BAM files for visualization in a genome browser like IGV.

2. Prerequisites

You must have the following software and R packages installed and available in your environment.
Software

    Snakemake: The pipeline management system.

    FastQC: For initial read quality assessment.

    MultiQC: To aggregate FastQC reports.

    Trimmomatic: For adapter and quality trimming.

    Trinity: For de novo transcriptome assembly.

    Kallisto: For fast read quantification.

    Bowtie2: For aligning reads to the Trinity assembly.

    Samtools: For sorting and indexing BAM files.

    R: The statistical computing environment.

    A PBS-based HPC cluster (optional): For running the pipeline on a cluster.

R Packages

The scripts/rna_analysis.R script requires the following R packages, which you should install before running the analysis step.

install.packages(c("tidyverse", "ggplot2", "pheatmap", "plotly", "yaml"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("tximport", "DESeq2", "PCAtools", "GenomicRanges"))

3. Directory Structure

The pipeline expects the following directory structure:

.
├── Snakefile
├── config.yaml
├── submit.sh
├── scripts
│   └── rna_analysis.R
└── data
    └── raw
        ├── ctrl_rep1_R1.fastq.gz
        └── ...

The pipeline will create a results/ directory to store all outputs:

.
└── results
    ├── fastqc_raw
    ├── multiqc
    ├── trimmed
    ├── trinity
    ├── bowtie2_index
    ├── bowtie2_align
    ├── kallisto_index
    ├── counts
    └── R_analysis

4. Configuration

The config.yaml file is the central place to configure the pipeline. You must edit this file to match your data and system resources.

    samples: Define your samples and the paths to their paired-end FASTQ files.

    threads: The number of CPU threads to use for each rule.

    trinity_memory: The maximum memory to allocate to Trinity (e.g., "30G").

    use_scratchdir: Set to true when running on an HPC to enable fast I/O using the temporary scratch directory ($SCRATCHDIR). Set to false for local runs.

Example config.yaml:

samples:
    ctrl_rep1:
        R1: "data/raw/ctrl_rep1_R1.fastq.gz"
        R2: "data/raw/ctrl_rep1_R2.fastq.gz"
    # ... more samples
threads: 8
trinity_memory: "30G"
use_scratchdir: false

5. Usage
5.1. Local Run

If you are running the pipeline on a single machine, you can simply call Snakemake from the command line.

    Set use_scratchdir: false in your config.yaml.

    Run the pipeline with the desired number of cores:

    snakemake --cores 8

5.2. HPC Run (PBS)

For an HPC environment with a PBS scheduler, use the provided submit.sh script. This script will submit the Snakemake job to the cluster, which will then handle the submission of individual rules as sub-jobs.

    Ensure use_scratchdir: true in your config.yaml.

    Submit the job using qsub:

    qsub submit.sh

    The job will run in the background, and the output will be written to snakemake_output.log and snakemake_error.log.

6. Key Output Files

    results/multiqc/multiqc_report.html: An interactive report summarizing the QC metrics of all samples.

    results/trinity/Trinity.fasta: The final de novo assembled transcriptome.

    results/bowtie2_align/{sample}.bam.bai: The indexed BAM files for visualization.

    results/counts/{sample}_kallisto.tsv: The raw counts for each contig from Kallisto.

    results/R_analysis/: This directory will contain all the plots and the differential expression results table (deseq2_results.csv).

7. Customization

The scripts/rna_analysis.R script is designed to be easily modified. You can open this file to:

    Adjust parameters for DESeq2.

    Change the threshold for significance in the Volcano plot.

    Add new plots or statistical analyses (e.g., GSEA, GO enrichment).

    Change the color schemes or plot aesthetics.
