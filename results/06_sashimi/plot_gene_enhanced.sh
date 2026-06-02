#!/bin/bash
# Helper script to generate enhanced sashimi plot for a specific gene
# Usage: bash plot_gene_enhanced.sh GENE_NAME [REGION]
# Example: bash plot_gene_enhanced.sh Rbfox1
# Example: bash plot_gene_enhanced.sh Myc chr15:61985000-61992000

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 GENE_NAME [chr:start-end]"
    exit 1
fi

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/ggsashimi

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
OUTPUT_DIR="${BASE_DIR}/results/06_sashimi"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"
BAM_CONFIG="${OUTPUT_DIR}/bam_config_enhanced.tsv"

GENE=$1
REGION=$2

mkdir -p ${OUTPUT_DIR}/user_genes/${GENE}

if [[ -z "${REGION}" ]]; then
    REGION=$(grep -P "\tgene\t.*gene_name \"${GENE}\"" ${GTF_FILE} | head -1 | awk '{print $1":"$4"-"$5}')
    if [[ -z "${REGION}" ]]; then
        echo "Gene ${GENE} not found in GTF"
        exit 1
    fi
fi

echo "Generating enhanced sashimi plot for ${GENE} at ${REGION}"

# Expand region for context
CHR=$(echo ${REGION} | cut -d':' -f1)
START=$(echo ${REGION} | cut -d':' -f2 | cut -d'-' -f1)
END=$(echo ${REGION} | cut -d'-' -f2)
SPAN=$((END - START))
EXPAND=$((SPAN / 10))
NEW_START=$((START - EXPAND))
NEW_END=$((END + EXPAND))
if [[ ${NEW_START} -lt 1 ]]; then NEW_START=1; fi
EXPANDED_REGION="${CHR}:${NEW_START}-${NEW_END}"

# Generate PDF
ggsashimi.py \
    -b ${BAM_CONFIG} \
    -c ${EXPANDED_REGION} \
    -g ${GTF_FILE} \
    -o ${OUTPUT_DIR}/user_genes/${GENE}/${GENE}_sashimi \
    -M 50 \
    -C 3 \
    --alpha 0.25 \
    --height 5 \
    --width 18 \
    --base-size 12 \
    --ann-height 2 \
    --shrink \
    -F pdf

# Generate PNG
ggsashimi.py \
    -b ${BAM_CONFIG} \
    -c ${EXPANDED_REGION} \
    -g ${GTF_FILE} \
    -o ${OUTPUT_DIR}/user_genes/${GENE}/${GENE}_sashimi \
    -M 50 \
    -C 3 \
    --alpha 0.25 \
    --height 5 \
    --width 18 \
    --base-size 12 \
    --ann-height 2 \
    --shrink \
    -F png

echo "Done! Check ${OUTPUT_DIR}/user_genes/${GENE}/"
