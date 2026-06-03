#!/bin/bash
#===============================================================================
# 13_AS_genelists.sh  --  Foundation step for the AS interpretation sub-pipeline
#
# Parses the rMATS output for the two MUTATION comparisons (EMX1_wt_vs_mut,
# Nestin_wt_vs_mut), applies the n=1-defensible filter (|dPSI|>=0.1 + read
# support; NOT p-value), and emits the gene lists, the EMX1 n Nestin concordant
# set, and the rMATS-tested background that steps 14-17 consume.
#
# Run FIRST, before 14/15/16/17.
#
# USAGE:  sbatch 13_AS_genelists.sh
#===============================================================================
#SBATCH --job-name=13_AS_genelists
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/13_AS_genelists.err"
#SBATCH --output="./logs/13_AS_genelists.out"

set -euo pipefail

# Conda: reuse the scRNA analysis env (has pandas/numpy). Adjust if needed.
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate rna_seq_analysis_deep

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
SCRIPT_DIR="${BASE_DIR}/scripts"

mkdir -p "${BASE_DIR}/logs"

echo "=== Step 13: AS gene lists & concordance ==="
echo "Timestamp: $(date)"

python "${SCRIPT_DIR}/13_AS_genelists.py" \
    --base-dir "${BASE_DIR}" \
    --min-dpsi 0.10 \
    --min-reads 20

echo "=== Step 13 complete ==="
echo "Outputs: ${BASE_DIR}/results/09_AS_interpretation/genelists/"
echo "Timestamp: $(date)"
