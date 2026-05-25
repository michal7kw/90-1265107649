#!/bin/bash

#===============================================================================
# SCRIPT: 0_build_star_index.sh
# PURPOSE: Build STAR genome index for mouse GRCm39
#
# DESCRIPTION:
# Builds a STAR index compatible with STAR 2.7.11b for the mouse genome.
# This is required because the existing Cell Ranger index was built with
# STAR 2.7.1a which is incompatible.
#
# USAGE:
# sbatch 0_build_star_index.sh
#===============================================================================

#SBATCH --job-name=star_index_mouse
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/0_star_index.err"
#SBATCH --output="./logs/0_star_index.out"

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/rnaseq-quant

# Paths
GENOME_FASTA="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/fasta/genome.fa"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/STAR_GRCm39"

echo "=== Building STAR index for mouse GRCm39 ==="
echo "Timestamp: $(date)"
echo "STAR version: $(STAR --version)"
echo ""
echo "Genome FASTA: ${GENOME_FASTA}"
echo "GTF file: ${GTF_FILE}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Build STAR index
# Using sjdbOverhang=100 (read length - 1, for 150bp reads this is fine)
STAR \
    --runMode genomeGenerate \
    --runThreadN 16 \
    --genomeDir ${OUTPUT_DIR} \
    --genomeFastaFiles ${GENOME_FASTA} \
    --sjdbGTFfile ${GTF_FILE} \
    --sjdbOverhang 100

echo ""
echo "=== STAR index generation complete ==="
echo "Timestamp: $(date)"
echo ""
echo "Index files:"
ls -lh ${OUTPUT_DIR}
