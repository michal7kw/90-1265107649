#!/bin/bash

#===============================================================================
# SCRIPT: 1_fastqc.sh
# PURPOSE: Quality Control for RiboTag RNA-seq Data
#
# DESCRIPTION:
# Performs quality control analysis of raw paired-end FASTQ files using FastQC.
# First step in the RiboTag analysis pipeline.
#
# USAGE:
# sbatch 1_fastqc.sh
#===============================================================================

#SBATCH --job-name=1_fastqc_ribotag
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --array=0-3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/1_fastqc_%a.err"
#SBATCH --output="./logs/1_fastqc_%a.out"

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/quality

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
RAW_DIR="${BASE_DIR}/00_fastq"
OUTPUT_DIR="${BASE_DIR}/results/01_fastqc"

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "=== Starting FastQC for sample: ${SAMPLE} ==="
echo "Timestamp: $(date)"
echo "Input directory: ${RAW_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Check if input files exist
R1_FILE="${RAW_DIR}/${SAMPLE}_R1_001.fastq.gz"
R2_FILE="${RAW_DIR}/${SAMPLE}_R2_001.fastq.gz"

if [[ ! -f "${R1_FILE}" ]]; then
    echo "ERROR: R1 file not found: ${R1_FILE}"
    exit 1
fi
if [[ ! -f "${R2_FILE}" ]]; then
    echo "ERROR: R2 file not found: ${R2_FILE}"
    exit 1
fi

echo "Input files found:"
echo "  R1: ${R1_FILE}"
echo "  R2: ${R2_FILE}"

echo "=== Running FastQC ==="
echo "Timestamp: $(date)"

# Run FastQC on both R1 and R2 files
fastqc \
    --threads 8 \
    --outdir ${OUTPUT_DIR} \
    --format fastq \
    ${R1_FILE} ${R2_FILE}

echo "=== Checking output files ==="
if [[ -f "${OUTPUT_DIR}/${SAMPLE}_R1_001_fastqc.html" && -f "${OUTPUT_DIR}/${SAMPLE}_R2_001_fastqc.html" ]]; then
    echo "SUCCESS: FastQC reports generated for ${SAMPLE}"
    echo "  R1 report: ${OUTPUT_DIR}/${SAMPLE}_R1_001_fastqc.html"
    echo "  R2 report: ${OUTPUT_DIR}/${SAMPLE}_R2_001_fastqc.html"
else
    echo "ERROR: FastQC reports not generated!"
    exit 1
fi

echo "=== FastQC complete for ${SAMPLE} ==="
echo "Timestamp: $(date)"
