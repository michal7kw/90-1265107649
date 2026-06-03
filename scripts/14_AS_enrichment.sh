#!/bin/bash
#===============================================================================
# 14_AS_enrichment.sh  --  GO/pathway over-representation on AS gene lists
#
# ORA (Enrichr via gseapy) on the step-13 gene lists, using the rMATS-tested
# gene universe as background and adding an RNA-processing/spliceosome focus
# pass. Prioritise the `concordant` set. Results are EXPLORATORY (n=1).
#
# PREREQUISITE: 13_AS_genelists.sh must have run.
# NOTE: needs internet access (Enrichr API) from the node, and `gseapy`.
#
# USAGE:  sbatch 14_AS_enrichment.sh
#===============================================================================
#SBATCH --job-name=14_AS_enrichment
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/14_AS_enrichment.err"
#SBATCH --output="./logs/14_AS_enrichment.out"

set -euo pipefail

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate rna_seq_analysis_deep   # has gseapy (used by 4_GO_Terms)

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"

mkdir -p "${BASE_DIR}/logs"

echo "=== Step 14: AS GO/pathway enrichment ==="
echo "Timestamp: $(date)"

# If gseapy is missing from the env, uncomment:
# pip install --quiet gseapy

python "${BASE_DIR}/scripts/14_AS_enrichment.py" --base-dir "${BASE_DIR}"

echo "=== Step 14 complete ==="
echo "Outputs: ${BASE_DIR}/results/09_AS_interpretation/enrichment/"
echo "Timestamp: $(date)"
