#!/bin/bash
# Helper script to generate sashimi plot for a specific gene
# Usage: bash plot_gene.sh GENE_NAME [REGION]
# Example: bash plot_gene.sh TP53
# Example: bash plot_gene.sh MYC chr8:127735434-127742951

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 GENE_NAME [chr:start-end]"
    exit 1
fi

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/ggsashimi

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
OUTPUT_DIR="${BASE_DIR}/results/06_sashimi"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"
BAM_CONFIG="${OUTPUT_DIR}/bam_config.tsv"

GENE=$1
REGION=$2

mkdir -p ${OUTPUT_DIR}/${GENE}

if [[ -z "${REGION}" ]]; then
    REGION=$(grep -P "\tgene\t.*gene_name \"${GENE}\"" ${GTF_FILE} | head -1 | awk '{print $1":"$4"-"$5}')
    if [[ -z "${REGION}" ]]; then
        echo "Gene ${GENE} not found in GTF"
        exit 1
    fi
fi

echo "Generating sashimi plot for ${GENE} at ${REGION}"

ggsashimi.py \
    -b ${BAM_CONFIG} \
    -c ${REGION} \
    -g ${GTF_FILE} \
    -o ${OUTPUT_DIR}/${GENE}/${GENE}_sashimi \
    -M 10 \
    -C 3 \
    --alpha 0.25 \
    --height 2 \
    --width 8 \
    --shrink \
    -F pdf

ggsashimi.py \
    -b ${BAM_CONFIG} \
    -c ${REGION} \
    -g ${GTF_FILE} \
    -o ${OUTPUT_DIR}/${GENE}/${GENE}_sashimi \
    -M 10 \
    -C 3 \
    --alpha 0.25 \
    --height 2 \
    --width 8 \
    --shrink \
    -F png

echo "Done! Check ${OUTPUT_DIR}/${GENE}/"
