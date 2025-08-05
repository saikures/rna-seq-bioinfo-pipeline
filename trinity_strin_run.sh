#!/bin/bash
#PBS -N trinity_assembly
#PBS -M saikures@gmail.com
#PBS -m abe
#PBS -l select=1:ncpus=32:mem=1500gb:scratch_local=1500gb
#PBS -l walltime=48:00:00
#PBS -j oe

# Load conda
source /storage/brno12-cerit/home/saikures/miniconda3/bin/activate
conda activate trinity_env

#trap "clean_scratch" TERM EXIT

# Define paths
REMOTE_USER="saikures"
REMOTE_HOST="storage-brno12-cerit.metacentrum.cz"
REMOTE_DATA_PATH="~/saikures_main/RNASeq_major/maj_dataset/trim_maj_stringent/"
FINALOUTDIR="$REMOTE_DATA_PATH/scratch_trinity_maj_strin_final"

mkdir -p "$FINALOUTDIR"

# Check scratch
if [ -z "${SCRATCHDIR:-}" ]; then
    echo "❌ SCRATCHDIR is not set!"
    exit 1
fi

cd "$SCRATCHDIR"
echo "$PBS_JOBID running on $(hostname -f) using $SCRATCHDIR" >> "$FINALOUTDIR/jobs_info.txt"

# Prepare read lists
LEFT_LIST=""
RIGHT_LIST=""

export TMPDIR=$SCRATCHDIR

shopt -s nullglob
for LEFT_REMOTE in $(ssh ${REMOTE_USER}@${REMOTE_HOST} "ls ${REMOTE_DATA_PATH}/M-*_1_paired.fastq.gz"); do
    SAMPLE_ID=$(basename "${LEFT_REMOTE}" _1_paired.fastq.gz)
    RIGHT_REMOTE="${REMOTE_DATA_PATH}/${SAMPLE_ID}_2_paired.fastq.gz"

    # Check right read exists remotely
    if ssh ${REMOTE_USER}@${REMOTE_HOST} "[ -f ${RIGHT_REMOTE} ]"; then
        echo "✅ $SAMPLE_ID"
        
        # Copy both files to SCRATCHDIR
        scp "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_DATA_PATH}/${SAMPLE_ID}_1_paired.fastq.gz" .
        scp "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_DATA_PATH}/${SAMPLE_ID}_2_paired.fastq.gz" .

        LEFT_LIST+="$SCRATCHDIR/${SAMPLE_ID}_1_paired.fastq.gz,"
        RIGHT_LIST+="$SCRATCHDIR/${SAMPLE_ID}_2_paired.fastq.gz,"
    else
        echo "⚠️  Skipping $SAMPLE_ID: right read not found"
    fi
done
shopt -u nullglob

# Remove trailing commas
LEFT_LIST=${LEFT_LIST%,}
RIGHT_LIST=${RIGHT_LIST%,}

echo "📄 Files in scratch:"
ls -lh *.fastq.gz

# Run Trinity
echo "🚀 Running Trinity..."
Trinity \
  --seqType fq \
  --max_memory 1500G \
  --min_contig_length 200 \
  --CPU 32 \
  --left "$LEFT_LIST" \
  --right "$RIGHT_LIST" \
  --output "$SCRATCHDIR/trinity_out" \
  > trinity_run.log 2> trinity_run.err

scp -r "$SCRATCHDIR/trinity_out/"* "${REMOTE_USER}@${REMOTE_HOST}:$FINALOUTDIR/" || export CLEAN_SCRATCH=false


