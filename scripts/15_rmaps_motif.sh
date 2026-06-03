#!/bin/bash
#===============================================================================
# 15_rmaps_motif.sh  --  RBP / splicing-factor motif maps around regulated exons
#
# Identifies WHICH RNA-binding proteins may drive the observed splicing changes
# by mapping RBP motif enrichment around the regulated exons. This is the
# "splicing-factor regulatory network" step: because Mbnl1 / Hnrnpa2b1 are
# themselves AS hits, enrichment of their motifs here tests the cascade /
# feedback hypothesis directly.
#
# PRIMARY TOOL: rMAPS2 (standalone) -- purpose-built to take rMATS MATS output
#   and emit per-nucleotide RBP motif density maps for SE/A3SS/A5SS.
#   Download: http://rmaps.cecsresearch.org/  (standalone package)
#   Point RMAPS2_DIR at the install (containing the rMAPS2 launcher).
#
# FALLBACK (always runnable): 15_rbp_motif_fallback.py -- known-motif Fisher
#   enrichment in intronic windows flanking SE exons (fg=regulated vs bg=tested).
#   Used automatically when RMAPS2_DIR is unset or the launcher is missing.
#
# PREREQUISITE: 13_AS_genelists.sh (for the filter parameters / context).
# Operates on the raw rMATS files in results/05_splicing/<comp>/.
#
# USAGE:  sbatch 15_rmaps_motif.sh
#===============================================================================
#SBATCH --job-name=15_rmaps_motif
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/15_rmaps_motif.err"
#SBATCH --output="./logs/15_rmaps_motif.out"

set -euo pipefail

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate rna_seq_analysis_deep   # needs pyfaidx, scipy for the fallback

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
SPLICING_DIR="${BASE_DIR}/results/05_splicing"
OUT_DIR="${BASE_DIR}/results/09_AS_interpretation/rmaps"
GENOME_FA="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/fasta/genome.fa"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"

# Set this if you have installed the rMAPS2 standalone:
RMAPS2_DIR="${RMAPS2_DIR:-}"

mkdir -p "${OUT_DIR}" "${BASE_DIR}/logs"

echo "=== Step 15: RBP motif maps ==="
echo "Timestamp: $(date)"

COMPARISONS=("EMX1_wt_vs_mut" "Nestin_wt_vs_mut")

if [[ -n "${RMAPS2_DIR}" && -e "${RMAPS2_DIR}/rmaps2.py" ]]; then
    echo ">> rMAPS2 found at ${RMAPS2_DIR}; running motif maps."
    for COMP in "${COMPARISONS[@]}"; do
        for EVENT in SE A3SS A5SS; do
            MATS="${SPLICING_DIR}/${COMP}/${EVENT}.MATS.JC.txt"
            COMP_OUT="${OUT_DIR}/${COMP}/${EVENT}"
            mkdir -p "${COMP_OUT}"
            echo "   rMAPS2: ${COMP} ${EVENT}"
            # NOTE: confirm the exact flag names against your rMAPS2 version;
            # the standalone expects the rMATS MATS file, an event type, the
            # genome FASTA, and a significance definition for fg vs bg.
            python "${RMAPS2_DIR}/rmaps2.py" \
                --splicing "${MATS}" \
                --event "${EVENT}" \
                --fasta "${GENOME_FA}" \
                --gtf "${GTF_FILE}" \
                --dpsi 0.1 --fdr 1.0 \
                --output "${COMP_OUT}" \
                || echo "   [warn] rMAPS2 failed for ${COMP}/${EVENT}; see log."
        done
    done
else
    echo ">> rMAPS2 not configured (RMAPS2_DIR unset/missing)."
    echo ">> Running self-contained fallback motif enrichment instead."
    if [[ ! -e "${GENOME_FA}" ]]; then
        echo "ERROR: genome FASTA not found at ${GENOME_FA}" >&2
        echo "       Set GENOME_FA to the GRCm39 primary assembly FASTA." >&2
        exit 1
    fi
    # pyfaidx is required by the fallback:
    python -c "import pyfaidx" 2>/dev/null || pip install --quiet pyfaidx
    python "${BASE_DIR}/scripts/15_rbp_motif_fallback.py" \
        --base-dir "${BASE_DIR}" \
        --genome "${GENOME_FA}"
fi

echo "=== Step 15 complete ==="
echo "Outputs: ${OUT_DIR}/"
echo "Timestamp: $(date)"
