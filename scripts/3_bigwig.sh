#!/bin/bash

#===============================================================================
# SCRIPT: 3_bigwig.sh (CORRECTED VERSION)
# PURPOSE: Generate BigWig Coverage Tracks for RiboTag RNA-seq Data
#
# DESCRIPTION:
# Creates normalized BigWig coverage tracks from aligned BAM files using
# deepTools bamCoverage. These tracks can be visualized in genome browsers
# like IGV, UCSC Genome Browser, or JBrowse.
#
# CORRECTIONS APPLIED:
# 1. Removed --extendReads flag (was causing signal across splice junctions)
# 2. Added --minMappingQuality 255 (filters to uniquely mapped reads only)
#
# OUTPUT:
# - CPM-normalized BigWig files (for comparing between samples)
# - Forward and reverse strand-specific tracks (optional)
#
# USAGE:
# cd bigwig_corrected && sbatch scripts/3_bigwig.sh
#===============================================================================

#SBATCH --job-name=3_bigwig_corrected
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=0-3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/3_bigwig_%a.err"
#SBATCH --output="./logs/3_bigwig_%a.out"

# Set up conda environment with deeptools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/bigwig

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
ALIGNED_DIR="${BASE_DIR}/results/02_aligned"
OUTPUT_DIR="${BASE_DIR}/bigwig_corrected/results/03_bigwig"

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Sample name mapping for better output file names
declare -A SAMPLE_NAMES
SAMPLE_NAMES[1]="EMX1_hippo_wt"
SAMPLE_NAMES[2]="EMX1_hippo_mut"
SAMPLE_NAMES[3]="Nestin_hippo_wt"
SAMPLE_NAMES[4]="Nestin_hippo_mut"

SAMPLE_NAME=${SAMPLE_NAMES[$SAMPLE]}

echo "=== Starting CORRECTED BigWig generation for sample: ${SAMPLE} (${SAMPLE_NAME}) ==="
echo "Timestamp: $(date)"
echo ""
echo "CORRECTIONS APPLIED:"
echo "  - Removed --extendReads (prevents signal across splice junctions)"
echo "  - Added --minMappingQuality 255 (uniquely mapped reads only)"
echo ""

# Input BAM file
BAM_FILE="${ALIGNED_DIR}/${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam"

if [[ ! -f "${BAM_FILE}" ]]; then
    echo "ERROR: BAM file not found: ${BAM_FILE}"
    exit 1
fi

if [[ ! -f "${BAM_FILE}.bai" ]]; then
    echo "ERROR: BAM index not found. Indexing now..."
    samtools index -@ 16 ${BAM_FILE}
fi

echo "Input BAM: ${BAM_FILE}"
echo "Output directory: ${OUTPUT_DIR}"

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "=== Step 1: Generating CPM-normalized BigWig (unstranded) ==="
echo "Timestamp: $(date)"

# Generate CPM-normalized BigWig (counts per million)
# This is the main track for visualization
bamCoverage \
    --bam ${BAM_FILE} \
    --outFileName ${OUTPUT_DIR}/${SAMPLE_NAME}_CPM.bw \
    --outFileFormat bigwig \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 16 \
    --minMappingQuality 255 \
    --ignoreDuplicates

if [[ -f "${OUTPUT_DIR}/${SAMPLE_NAME}_CPM.bw" ]]; then
    echo "SUCCESS: CPM-normalized BigWig created"
    ls -lh ${OUTPUT_DIR}/${SAMPLE_NAME}_CPM.bw
else
    echo "ERROR: Failed to create CPM-normalized BigWig"
    exit 1
fi

echo "=== Step 2: Generating strand-specific BigWig tracks ==="
echo "Timestamp: $(date)"

# Forward strand (for genes on + strand)
bamCoverage \
    --bam ${BAM_FILE} \
    --outFileName ${OUTPUT_DIR}/${SAMPLE_NAME}_forward.bw \
    --outFileFormat bigwig \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 16 \
    --filterRNAstrand forward \
    --minMappingQuality 255 \
    --ignoreDuplicates

# Reverse strand (for genes on - strand)
bamCoverage \
    --bam ${BAM_FILE} \
    --outFileName ${OUTPUT_DIR}/${SAMPLE_NAME}_reverse.bw \
    --outFileFormat bigwig \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 16 \
    --filterRNAstrand reverse \
    --minMappingQuality 255 \
    --ignoreDuplicates

echo "=== Step 3: Generating RPKM-normalized BigWig ==="
echo "Timestamp: $(date)"

# RPKM normalization (alternative normalization)
bamCoverage \
    --bam ${BAM_FILE} \
    --outFileName ${OUTPUT_DIR}/${SAMPLE_NAME}_RPKM.bw \
    --outFileFormat bigwig \
    --normalizeUsing RPKM \
    --binSize 10 \
    --numberOfProcessors 16 \
    --minMappingQuality 255 \
    --ignoreDuplicates

echo "=== Checking output files ==="
echo "Generated BigWig files for ${SAMPLE_NAME}:"
ls -lh ${OUTPUT_DIR}/${SAMPLE_NAME}*.bw

echo "=== CORRECTED BigWig generation complete for ${SAMPLE} ==="
echo "Timestamp: $(date)"
