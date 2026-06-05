#!/bin/bash
#===============================================================================
# 17_AS_integration.sh  --  External integration + interpretation report
#
# Overlaps the high-confidence concordant AS genes (step 13) with the snRNA SJU
# convergence table, the Soelter et al. 2026 RBPs, curated disease gene sets, and
# a STRING PPI network, then writes the interpretation report
# (results/09_AS_interpretation/AS_interpretation_report.md). Folds in the
# overview files from steps 14/15/16 if they exist.
#
# PREREQUISITE: 13_AS_genelists.sh (required); 14/15/16 optional but recommended.
# NOTE: STRING needs internet; pass --no-string to skip.
#
# USAGE:  sbatch 17_AS_integration.sh
#===============================================================================
#SBATCH --job-name=17_AS_integration
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/17_AS_integration.err"
#SBATCH --output="./logs/17_AS_integration.out"

set -euo pipefail

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate rna_seq_analysis_deep   # needs scipy; requests for STRING (opt.)

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
RNA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA"

mkdir -p "${BASE_DIR}/logs"

echo "=== Step 17: AS external integration + report ==="
echo "Timestamp: $(date)"

python "${BASE_DIR}/scripts/17_AS_integration.py" \
    --base-dir "${BASE_DIR}" \
    --rna-dir "${RNA_DIR}"

echo "=== Step 17 complete ==="
echo "Report: ${BASE_DIR}/results/09_AS_interpretation/AS_interpretation_report.md"
echo "Timestamp: $(date)"
