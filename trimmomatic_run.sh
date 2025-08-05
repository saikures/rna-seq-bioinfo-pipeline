#!/bin/bash
#PBS -N trim_maj_strin
#PBS -l select=1:ncpus=126:mem=256gb
#PBS -l walltime=24:00:00
#PBS -o trim_maj_strin.log
#PBS -j oe

export OMP_NUM_THREADS=126

TRIM_JAR=/storage/brno2/home/saikures/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/trimmomatic.jar
ADAPTERS=/storage/brno2/home/saikures/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa
THREADS=126
OUTDIR=/storage/brno12-cerit/home/saikures/saikures_main/RNASeq_major/maj_dataset/trim_maj_stringent
FASTP_TMPDIR="$OUTDIR/fastp_polyG_trimmed"

cd "$PBS_O_WORKDIR"

mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR/fastp_polyG_trimmed"

module load fastp
module load trimmomatic
module load jdk

for R1 in M-*_1.fastq.gz; do
    SAMPLE_BASE=$(basename "$R1" _1.fastq.gz)
    R2="${SAMPLE_BASE}_2.fastq.gz"

    # Final trimmed paired output file to check
    FINAL_TRIMMED_FILE="$OUTDIR/${SAMPLE_BASE}_1_paired.fastq.gz"

    # Skip if final output file already exists
    if [[ -f "$FINAL_TRIMMED_FILE" ]]; then
        echo "Skipping $SAMPLE_BASE — already trimmed."
        continue
    fi

    if [[ -f "$R2" ]]; then
        echo "Processing paired-end sample: $SAMPLE_BASE"

        # Step 1: Poly-G trimming with fastp
        fastp \
            -i "$R1" -I "$R2" \
            -o "$FASTP_TMPDIR/${SAMPLE_BASE}_1.fastq.gz" \
            -O "$FASTP_TMPDIR/${SAMPLE_BASE}_2.fastq.gz" \
            --trim_poly_g \
            --thread $THREADS \

        # Step 2: Adapter and quality trimming with Trimmomatic
        java -jar "$TRIM_JAR" PE -threads $THREADS -phred33 \
            "$FASTP_TMPDIR/${SAMPLE_BASE}_1.fastq.gz" "$FASTP_TMPDIR/${SAMPLE_BASE}_2.fastq.gz" \
            "$OUTDIR/${SAMPLE_BASE}_1_paired.fastq.gz" "$OUTDIR/${SAMPLE_BASE}_1_unpaired.fastq.gz" \
            "$OUTDIR/${SAMPLE_BASE}_2_paired.fastq.gz" "$OUTDIR/${SAMPLE_BASE}_2_unpaired.fastq.gz" \
            ILLUMINACLIP:"$ADAPTERS":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
    else
        echo "WARNING: Missing $R2 for $SAMPLE_BASE — skipping..."
    fi
done
