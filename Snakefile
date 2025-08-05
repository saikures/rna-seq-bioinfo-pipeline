# ==============================================================================
# SNAKEMAKE CONFIGURATION
# ==============================================================================
configfile: "config.yaml"

# ==============================================================================
# WILDCARDS
# ==============================================================================
SAMPLES = list(config["samples"].keys())

# ==============================================================================
# PIPELINE RULES
# ==============================================================================

# ------------------------------------------------------------------------------
# Rule 0: All
# This is the main rule that triggers all other rules.
# It defines the final output files of the entire pipeline.
# ------------------------------------------------------------------------------
rule all:
    input:
        "results/multiqc/multiqc_report.html",
        "results/trinity/Trinity.fasta",
        expand("results/bowtie2_align/{sample}.bam.bai", sample=SAMPLES),
        expand("results/counts/{sample}_kallisto.tsv", sample=SAMPLES),
        "results/R_analysis/deseq2_results.csv",
        "results/R_analysis/pca_plot.png",
        "results/R_analysis/ma_plot.png",
        "results/R_analysis/volcano_plot.png",
        "results/R_analysis/p_value_plot.png",
        "results/R_analysis/gene_regulation_plot.png",
        "results/R_analysis/heatmap.png"


# ------------------------------------------------------------------------------
# Rule 1: Quality Control (QC)
# ------------------------------------------------------------------------------
# QC on raw data (no scratchdir needed here as it's a quick analysis)
rule fastqc_raw:
    input:
        R1 = config["samples"]["{sample}"]["R1"],
        R2 = config["samples"]["{sample}"]["R2"]
    output:
        "results/fastqc_raw/{sample}_R1_fastqc.html",
        "results/fastqc_raw/{sample}_R2_fastqc.html"
    shell:
        "fastqc {input.R1} {input.R2} -o results/fastqc_raw/"

# multiqc on raw data
rule multiqc_raw:
    input:
        expand("results/fastqc_raw/{sample}_R1_fastqc.zip", sample=SAMPLES)
    output:
        "results/multiqc/multiqc_report_raw.html"
    shell:
        "multiqc results/fastqc_raw/ -o results/multiqc/"

# ------------------------------------------------------------------------------
# Rule 2: Trimming with Scratchdir Logic
# This rule copies files to scratch, runs Trimmomatic, and copies back.
# ------------------------------------------------------------------------------
rule trimmomatic:
    input:
        R1 = config["samples"]["{sample}"]["R1"],
        R2 = config["samples"]["{sample}"]["R2"]
    output:
        R1_paired = "results/trimmed/{sample}_R1.paired.fastq.gz",
        R1_unpaired = "results/trimmed/{sample}_R1.unpaired.fastq.gz",
        R2_paired = "results/trimmed/{sample}_R2.paired.fastq.gz",
        R2_unpaired = "results/trimmed/{sample}_R2.unpaired.fastq.gz"
    params:
        adapter_file = "TruSeq3-PE.fa" # You need to provide this file.
    threads: config["threads"]
    resources:
        # Define resources for the HPC scheduler, e.g., memory, time
        mem_mb = 4096
    shell:
        """
        if {config[use_scratchdir]}; then
            # Use a temporary directory on the scratch space
            SCRATCH_TEMP_DIR=$SCRATCHDIR/snakemake/trimmed/{wildcards.sample}
            mkdir -p $SCRATCH_TEMP_DIR
            echo "Copying raw data to scratch directory: $SCRATCH_TEMP_DIR"
            cp {input.R1} {input.R2} $SCRATCH_TEMP_DIR/

            # Run Trimmomatic using the files in the scratch directory
            trimmomatic PE -phred33 $SCRATCH_TEMP_DIR/{wildcards.sample}_R1.fastq.gz \
            $SCRATCH_TEMP_DIR/{wildcards.sample}_R2.fastq.gz \
            $SCRATCH_TEMP_DIR/{wildcards.sample}_R1.paired.fastq.gz \
            $SCRATCH_TEMP_DIR/{wildcards.sample}_R1.unpaired.fastq.gz \
            $SCRATCH_TEMP_DIR/{wildcards.sample}_R2.paired.fastq.gz \
            $SCRATCH_TEMP_DIR/{wildcards.sample}_R2.unpaired.fastq.gz \
            ILLUMINACLIP:{params.adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

            # Copy the results back to the permanent storage
            echo "Copying trimmed results back to permanent storage."
            cp $SCRATCH_TEMP_DIR/{wildcards.sample}_R1.paired.fastq.gz {output.R1_paired}
            cp $SCRATCH_TEMP_DIR/{wildcards.sample}_R2.paired.fastq.gz {output.R2_paired}
            cp $SCRATCH_TEMP_DIR/{wildcards.sample}_R1.unpaired.fastq.gz {output.R1_unpaired}
            cp $SCRATCH_TEMP_DIR/{wildcards.sample}_R2.unpaired.fastq.gz {output.R2_unpaired}

            # Clean up the scratch directory
            echo "Cleaning up scratch directory."
            rm -r $SCRATCH_TEMP_DIR
        else
            # Run Trimmomatic directly on permanent storage
            trimmomatic PE -phred33 {input.R1} {input.R2} \
            {output.R1_paired} {output.R1_unpaired} \
            {output.R2_paired} {output.R2_unpaired} \
            ILLUMINACLIP:{params.adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        fi
        """

# ------------------------------------------------------------------------------
# Rule 3: De Novo Assembly with Trinity (Scratchdir is critical here)
# ------------------------------------------------------------------------------
rule trinity:
    input:
        expand("results/trimmed/{sample}_R1.paired.fastq.gz", sample=SAMPLES),
        expand("results/trimmed/{sample}_R2.paired.fastq.gz", sample=SAMPLES)
    output:
        "results/trinity/Trinity.fasta"
    params:
        memory = config["trinity_memory"]
    threads: config["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: int(config["trinity_memory"].replace("G", "")) * 1024
    shell:
        """
        if {config[use_scratchdir]}; then
            # Create a Trinity-specific scratch directory
            TRINITY_SCRATCH_DIR=$SCRATCHDIR/snakemake/trinity
            mkdir -p $TRINITY_SCRATCH_DIR

            # Copy trimmed reads to scratch for faster I/O
            echo "Copying trimmed reads to scratch directory for Trinity."
            cp {input} $TRINITY_SCRATCH_DIR/

            # Prepare left and right read lists for Trinity
            LEFT_READS=$(find $TRINITY_SCRATCH_DIR -name "*_R1.paired.fastq.gz" | tr '\\n' ',' | sed 's/,$//')
            RIGHT_READS=$(find $TRINITY_SCRATCH_DIR -name "*_R2.paired.fastq.gz" | tr '\\n' ',' | sed 's/,$//')

            # Run Trinity using the scratch directory
            Trinity --seqType fq --left $LEFT_READS --right $RIGHT_READS \
            --CPU {threads} --max_memory {params.memory} --output $TRINITY_SCRATCH_DIR/trinity_output

            # Copy the final assembly back to permanent storage
            echo "Copying final Trinity.fasta back to permanent storage."
            cp $TRINITY_SCRATCH_DIR/trinity_output/Trinity.fasta {output}

            # Clean up the scratch directory
            echo "Cleaning up scratch directory."
            rm -r $TRINITY_SCRATCH_DIR
        else
            LEFT_READS=$(echo {input} | cut -d' ' -f1,3)
            RIGHT_READS=$(echo {input} | cut -d' ' -f2,4)
            # Run Trinity directly on permanent storage
            Trinity --seqType fq --left $LEFT_READS --right $RIGHT_READS \
            --CPU {threads} --max_memory {params.memory} --output results/trinity/
        fi
        """

# ------------------------------------------------------------------------------
# Rule 4: Alignment to the Trinity Assembly (using scratchdir)
# ------------------------------------------------------------------------------
rule bowtie2_index:
    input:
        "results/trinity/Trinity.fasta"
    output:
        "results/bowtie2_index/Trinity.1.bt2"
    shell:
        "bowtie2-build {input} results/bowtie2_index/Trinity"

rule bowtie2_align:
    input:
        index_prefix = "results/bowtie2_index/Trinity",
        R1 = "results/trimmed/{sample}_R1.paired.fastq.gz",
        R2 = "results/trimmed/{sample}_R2.paired.fastq.gz"
    output:
        "results/bowtie2_align/{sample}.bam"
    params:
        threads = config["threads"]
    threads: config["threads"]
    resources:
        mem_mb = 8192
    shell:
        """
        if {config[use_scratchdir]}; then
            # Create a temporary directory on the scratch space
            SCRATCH_TEMP_DIR=$SCRATCHDIR/snakemake/bowtie2/{wildcards.sample}
            mkdir -p $SCRATCH_TEMP_DIR
            echo "Copying files to scratch directory: $SCRATCH_TEMP_DIR"
            cp {input.R1} {input.R2} $SCRATCH_TEMP_DIR/
            cp results/bowtie2_index/Trinity.1.bt2 $SCRATCH_TEMP_DIR/
            cp results/bowtie2_index/Trinity.2.bt2 $SCRATCH_TEMP_DIR/
            # ... repeat for all bowtie2 index files (3.bt2, 4.bt2, rev.1.bt2, rev.2.bt2)
            cp results/bowtie2_index/Trinity.3.bt2 $SCRATCH_TEMP_DIR/
            cp results/bowtie2_index/Trinity.4.bt2 $SCRATCH_TEMP_DIR/
            cp results/bowtie2_index/Trinity.rev.1.bt2 $SCRATCH_TEMP_DIR/
            cp results/bowtie2_index/Trinity.rev.2.bt2 $SCRATCH_TEMP_DIR/

            # Run bowtie2 using the scratch files
            bowtie2 -x $SCRATCH_TEMP_DIR/Trinity -1 $SCRATCH_TEMP_DIR/{wildcards.sample}_R1.paired.fastq.gz \
            -2 $SCRATCH_TEMP_DIR/{wildcards.sample}_R2.paired.fastq.gz -p {threads} | \
            samtools view -bS - | samtools sort -o $SCRATCH_TEMP_DIR/{wildcards.sample}.bam

            # Copy the results back to permanent storage
            echo "Copying aligned BAM back to permanent storage."
            cp $SCRATCH_TEMP_DIR/{wildcards.sample}.bam {output}
            
            # Clean up the scratch directory
            echo "Cleaning up scratch directory."
            rm -r $SCRATCH_TEMP_DIR
        else
            # Run bowtie2 directly on permanent storage
            bowtie2 -x {input.index_prefix} -1 {input.R1} -2 {input.R2} -p {threads} | \
            samtools view -bS - | samtools sort -o {output}
        fi
        """

# ------------------------------------------------------------------------------
# Rule 5: Quantification with Kallisto (using scratchdir)
# ------------------------------------------------------------------------------
rule kallisto_index:
    input:
        "results/trinity/Trinity.fasta"
    output:
        "results/kallisto_index/Trinity.idx"
    shell:
        "kallisto index -i {output} {input}"

rule kallisto_quant:
    input:
        index = "results/kallisto_index/Trinity.idx",
        R1 = "results/trimmed/{sample}_R1.paired.fastq.gz",
        R2 = "results/trimmed/{sample}_R2.paired.fastq.gz"
    output:
        "results/counts/{sample}_kallisto.tsv"
    params:
        threads = config["threads"]
    threads: config["threads"]
    resources:
        mem_mb = 8192
    shell:
        """
        if {config[use_scratchdir]}; then
            # Create a Kallisto-specific scratch directory
            KALLISTO_SCRATCH_DIR=$SCRATCHDIR/snakemake/kallisto/{wildcards.sample}
            mkdir -p $KALLISTO_SCRATCH_DIR
            
            # Copy necessary files
            echo "Copying files to scratch directory for Kallisto."
            cp {input.index} $KALLISTO_SCRATCH_DIR/
            cp {input.R1} $KALLISTO_SCRATCH_DIR/
            cp {input.R2} $KALLISTO_SCRATCH_DIR/
            
            # Run Kallisto using the scratch files
            kallisto quant -i $KALLISTO_SCRATCH_DIR/Trinity.idx -o $KALLISTO_SCRATCH_DIR/output_quant \
            -b 100 -t {threads} $KALLISTO_SCRATCH_DIR/{wildcards.sample}_R1.paired.fastq.gz \
            $KALLISTO_SCRATCH_DIR/{wildcards.sample}_R2.paired.fastq.gz
            
            # Copy the results back to permanent storage
            echo "Copying kallisto results back to permanent storage."
            mkdir -p results/counts/{wildcards.sample}
            cp $KALLISTO_SCRATCH_DIR/output_quant/abundance.tsv {output}

            # Clean up the scratch directory
            echo "Cleaning up scratch directory."
            rm -r $KALLISTO_SCRATCH_DIR
        else
            # Run Kallisto directly on permanent storage
            kallisto quant -i {input.index} -o results/counts/{wildcards.sample} \
            -b 100 -t {threads} {input.R1} {input.R2} && \
            mv results/counts/{wildcards.sample}/abundance.tsv {output}
        fi
        """

# ------------------------------------------------------------------------------
# Rule 6: R-based Data Analysis
# ------------------------------------------------------------------------------
# The R script will handle all complex statistical and plotting tasks.
# This rule does not need a scratchdir as it primarily reads from disk and writes back.
rule r_analysis:
    input:
        expand("results/counts/{sample}_kallisto.tsv", sample=SAMPLES),
        config_file = "config.yaml"
    output:
        "results/R_analysis/deseq2_results.csv",
        "results/R_analysis/pca_plot.png",
        "results/R_analysis/ma_plot.png",
        "results/R_analysis/volcano_plot.png",
        "results/R_analysis/p_value_plot.png",
        "results/R_analysis/gene_regulation_plot.png",
        "results/R_analysis/heatmap.png"
    script:
        "scripts/rna_analysis.R"

# ------------------------------------------------------------------------------
# Rule 7: Visualization of .bam
# ------------------------------------------------------------------------------
# This rule sorts and indexes the BAM files for visualization in tools like IGV.
rule visualize_bam:
    input:
        "results/bowtie2_align/{sample}.bam"
    output:
        "results/bowtie2_align/{sample}.bam.bai"
    shell:
        "samtools index {input}"

