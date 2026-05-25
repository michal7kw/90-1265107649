#!/bin/bash

#===============================================================================
# SCRIPT: 7_sashimi_plots.sh
# PURPOSE: Generate Sashimi Plots for Splicing Visualization
#
# DESCRIPTION:
# Creates sashimi plots to visualize splicing events using ggsashimi.
# Can be run for specific genes of interest.
#
# USAGE:
# sbatch 7_sashimi_plots.sh                    # Use default genes
# sbatch 7_sashimi_plots.sh GENE1 GENE2 ...   # Specify genes
#
# Or interactively:
# bash 7_sashimi_plots.sh TP53 BRCA1 MYC
#===============================================================================

#SBATCH --job-name=7_sashimi_ribotag
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/7_sashimi.err"
#SBATCH --output="./logs/7_sashimi.out"

# Load conda
source /opt/common/tools/ric.cosr/miniconda3/bin/activate

# Check if ggsashimi environment exists, if not create it
SASHIMI_ENV="/beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/ggsashimi"

if [[ ! -d "${SASHIMI_ENV}" ]]; then
    echo "Creating ggsashimi conda environment..."
    conda create -y -p ${SASHIMI_ENV} -c bioconda -c conda-forge \
        python=3.9 \
        samtools \
        r-base \
        r-ggplot2 \
        r-gridextra \
        r-data.table

    conda activate ${SASHIMI_ENV}

    # Install ggsashimi from GitHub
    pip install git+https://github.com/guigolab/ggsashimi.git
fi

conda activate ${SASHIMI_ENV}

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
ALIGNED_DIR="${BASE_DIR}/results/02_aligned"
OUTPUT_DIR="${BASE_DIR}/results/06_sashimi"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "=== Starting Sashimi Plot Generation ==="
echo "Timestamp: $(date)"

#-------------------------------------------------------------------------------
# Create BAM configuration file for ggsashimi
#-------------------------------------------------------------------------------
BAM_CONFIG="${OUTPUT_DIR}/bam_config.tsv"

cat > ${BAM_CONFIG} << 'EOF'
EMX1_wt	${ALIGNED_DIR}/1/1_Aligned.sortedByCoord.out.bam	1
EMX1_mut	${ALIGNED_DIR}/2/2_Aligned.sortedByCoord.out.bam	2
Nestin_wt	${ALIGNED_DIR}/3/3_Aligned.sortedByCoord.out.bam	3
Nestin_mut	${ALIGNED_DIR}/4/4_Aligned.sortedByCoord.out.bam	4
EOF

# Replace variables in config
sed -i "s|\${ALIGNED_DIR}|${ALIGNED_DIR}|g" ${BAM_CONFIG}

echo "BAM configuration file created: ${BAM_CONFIG}"
cat ${BAM_CONFIG}

#-------------------------------------------------------------------------------
# Function to generate sashimi plot for a gene
#-------------------------------------------------------------------------------
generate_sashimi() {
    local GENE=$1
    local REGION=$2  # Optional: specific region (chr:start-end)

    echo ""
    echo "=== Generating sashimi plot for: ${GENE} ==="

    local GENE_OUTPUT="${OUTPUT_DIR}/${GENE}"
    mkdir -p ${GENE_OUTPUT}

    # If region not provided, try to get it from GTF
    if [[ -z "${REGION}" ]]; then
        # Extract gene coordinates from GTF
        REGION=$(grep -P "\tgene\t.*gene_name \"${GENE}\"" ${GTF_FILE} | \
            head -1 | \
            awk '{print $1":"$4"-"$5}')

        if [[ -z "${REGION}" ]]; then
            echo "WARNING: Gene ${GENE} not found in GTF. Skipping..."
            return 1
        fi
        echo "  Region from GTF: ${REGION}"
    fi

    # Run ggsashimi
    ggsashimi.py \
        -b ${BAM_CONFIG} \
        -c ${REGION} \
        -g ${GTF_FILE} \
        -o ${GENE_OUTPUT}/${GENE}_sashimi \
        -M 10 \
        -C 3 \
        --alpha 0.25 \
        --height 2 \
        --width 8 \
        --base-size 12 \
        --ann-height 2 \
        --shrink \
        -F pdf

    # Also generate PNG
    ggsashimi.py \
        -b ${BAM_CONFIG} \
        -c ${REGION} \
        -g ${GTF_FILE} \
        -o ${GENE_OUTPUT}/${GENE}_sashimi \
        -M 10 \
        -C 3 \
        --alpha 0.25 \
        --height 2 \
        --width 8 \
        --base-size 12 \
        --ann-height 2 \
        --shrink \
        -F png

    if [[ -f "${GENE_OUTPUT}/${GENE}_sashimi.pdf" ]]; then
        echo "  SUCCESS: Sashimi plot created"
        ls -lh ${GENE_OUTPUT}/${GENE}_sashimi.*
    else
        echo "  WARNING: Sashimi plot may have failed"
    fi
}

#-------------------------------------------------------------------------------
# Process genes
#-------------------------------------------------------------------------------

# Default example genes (can be overridden by command line arguments)
# Note: Using mouse gene names (sentence case) for GRCm39 reference
DEFAULT_GENES=(
    "Actb"      # Housekeeping gene - good control
    "Gapdh"     # Housekeeping gene - good control
    "Rbfox1"    # RNA-binding protein, involved in splicing
    "Ptbp1"     # Polypyrimidine tract binding protein
)

# Use command line arguments if provided, otherwise use defaults
if [[ $# -gt 0 ]]; then
    GENES=("$@")
else
    GENES=("${DEFAULT_GENES[@]}")
fi

echo ""
echo "Genes to process: ${GENES[@]}"
echo ""

for GENE in "${GENES[@]}"; do
    generate_sashimi "${GENE}"
done

#-------------------------------------------------------------------------------
# Create helper script for custom genes
#-------------------------------------------------------------------------------
HELPER_SCRIPT="${OUTPUT_DIR}/plot_gene.sh"
cat > ${HELPER_SCRIPT} << 'HELPER_EOF'
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
HELPER_EOF

chmod +x ${HELPER_SCRIPT}

echo ""
echo "=== Sashimi Plot Generation Complete ==="
echo ""
echo "To generate plots for additional genes, use:"
echo "  bash ${HELPER_SCRIPT} GENE_NAME"
echo ""
echo "Example:"
echo "  bash ${HELPER_SCRIPT} TP53"
echo "  bash ${HELPER_SCRIPT} MYC chr8:127735434-127742951"
echo ""
echo "Results in: ${OUTPUT_DIR}/"
echo "Timestamp: $(date)"
