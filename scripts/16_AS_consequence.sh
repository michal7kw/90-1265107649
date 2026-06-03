#!/bin/bash
#===============================================================================
# 16_AS_consequence.sh  --  Functional consequence of each splicing event
#
# Annotates every filtered event with frame impact (segment length %% 3),
# CDS-vs-UTR location (overlap with GTF CDS), and an NMD prediction, then
# sanity-checks the NMD prediction against the DEG log2FoldChange. Turns
# "gene X is spliced" into a mechanistic statement about the protein.
#
# PREREQUISITE: rMATS output (results/05_splicing) + DEG (results/04_DEG).
# Uses Bioconductor rtracklayer + GenomicRanges (standard). `maser` optional.
#
# USAGE:  sbatch 16_AS_consequence.sh
#===============================================================================
#SBATCH --job-name=16_AS_consequence
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/16_AS_consequence.err"
#SBATCH --output="./logs/16_AS_consequence.out"

set -euo pipefail

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
# An R env with Bioconductor (rtracklayer, GenomicRanges, optparse).
conda activate seurat_full2

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"

mkdir -p "${BASE_DIR}/logs"

echo "=== Step 16: AS functional consequence ==="
echo "Timestamp: $(date)"

Rscript "${BASE_DIR}/scripts/16_AS_consequence.R" \
    --base-dir "${BASE_DIR}" \
    --gtf "${GTF_FILE}" \
    --min-dpsi 0.10 \
    --min-reads 20

echo "=== Step 16 complete ==="
echo "Outputs: ${BASE_DIR}/results/09_AS_interpretation/consequence/"
echo "Timestamp: $(date)"
