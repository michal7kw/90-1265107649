#!/bin/bash

#===============================================================================
# SCRIPT: 7_sashimi_plots_enhanced.sh
# PURPOSE: Generate Enhanced Sashimi Plots for Splicing Visualization
#
# DESCRIPTION:
# Creates high-quality sashimi plots using ggsashimi with improved settings:
#   - Larger plot dimensions (18x5 inches)
#   - Higher PNG resolution (300 DPI)
#   - Stricter junction filtering (-M 50)
#   - Automatic selection of top differential splicing events
#   - Both PDF and PNG output formats
#
# MONOPLICATE DATA NOTE:
# This project has n=1 per condition. Unlike the reference project, we cannot
# aggregate replicates. Each sample is displayed individually.
#
# USAGE:
# sbatch 7_sashimi_plots_enhanced.sh                    # Use auto-selected genes
# sbatch 7_sashimi_plots_enhanced.sh GENE1 GENE2 ...   # Specify genes
#
# Or interactively:
# bash 7_sashimi_plots_enhanced.sh Rbfox1 Ptbp1
#===============================================================================

#SBATCH --job-name=7_sashimi_enhanced
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/7_sashimi_enhanced.err"
#SBATCH --output="./logs/7_sashimi_enhanced.out"

# Load conda
source /opt/common/tools/ric.cosr/miniconda3/bin/activate

# Check if ggsashimi environment exists
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
    pip install git+https://github.com/guigolab/ggsashimi.git
fi

conda activate ${SASHIMI_ENV}

#===============================================================================
# Configuration
#===============================================================================

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
ALIGNED_DIR="${BASE_DIR}/results/02_aligned"
SPLICING_DIR="${BASE_DIR}/results/05_splicing"
OUTPUT_DIR="${BASE_DIR}/results/06_sashimi"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"

# Enhanced plot settings
MIN_JUNCTION_COV=50     # Stricter filtering - show only major isoforms
PLOT_HEIGHT=5
PLOT_WIDTH=18
PNG_RESOLUTION=300

echo "============================================"
echo "Enhanced Sashimi Plot Generation"
echo "Start time: $(date)"
echo "============================================"
echo ""
echo "Settings:"
echo "  - Min junction coverage: ${MIN_JUNCTION_COV}"
echo "  - Plot dimensions: ${PLOT_WIDTH}x${PLOT_HEIGHT}"
echo "  - PNG resolution: ${PNG_RESOLUTION} DPI"
echo "  - Output formats: PDF + PNG"
echo ""

# Create output directories
mkdir -p ${OUTPUT_DIR}/{top_splicing,user_genes}

#===============================================================================
# Create BAM Configuration File
#===============================================================================

BAM_CONFIG="${OUTPUT_DIR}/bam_config_enhanced.tsv"

# For monoplicate data, each sample gets its own track (no aggregation)
cat > ${BAM_CONFIG} << EOF
EMX1_wt	${ALIGNED_DIR}/1/1_Aligned.sortedByCoord.out.bam	EMX1_wt
EMX1_mut	${ALIGNED_DIR}/2/2_Aligned.sortedByCoord.out.bam	EMX1_mut
Nestin_wt	${ALIGNED_DIR}/3/3_Aligned.sortedByCoord.out.bam	Nestin_wt
Nestin_mut	${ALIGNED_DIR}/4/4_Aligned.sortedByCoord.out.bam	Nestin_mut
EOF

echo "BAM configuration file created: ${BAM_CONFIG}"
cat ${BAM_CONFIG}
echo ""

#===============================================================================
# Function: Generate Enhanced Sashimi Plot
#===============================================================================

generate_sashimi_enhanced() {
    local GENE=$1
    local REGION=$2  # Optional: specific region (chr:start-end)
    local CATEGORY=$3  # Category subdirectory (top_splicing, user_genes)

    echo ""
    echo "=== Generating enhanced sashimi plot for: ${GENE} ==="

    local GENE_OUTPUT="${OUTPUT_DIR}/${CATEGORY}/${GENE}"
    mkdir -p ${GENE_OUTPUT}

    # If region not provided, extract from GTF
    if [[ -z "${REGION}" ]]; then
        REGION=$(grep -P "\tgene\t.*gene_name \"${GENE}\"" ${GTF_FILE} | \
            head -1 | \
            awk '{print $1":"$4"-"$5}')

        if [[ -z "${REGION}" ]]; then
            echo "WARNING: Gene ${GENE} not found in GTF. Skipping..."
            return 1
        fi
        echo "  Region from GTF: ${REGION}"
    else
        echo "  Using provided region: ${REGION}"
    fi

    # Parse chromosome and expand region slightly for context
    local CHR=$(echo ${REGION} | cut -d':' -f1)
    local START=$(echo ${REGION} | cut -d':' -f2 | cut -d'-' -f1)
    local END=$(echo ${REGION} | cut -d'-' -f2)

    # Expand region by 10% on each side for context
    local SPAN=$((END - START))
    local EXPAND=$((SPAN / 10))
    local NEW_START=$((START - EXPAND))
    local NEW_END=$((END + EXPAND))
    if [[ ${NEW_START} -lt 1 ]]; then NEW_START=1; fi
    local EXPANDED_REGION="${CHR}:${NEW_START}-${NEW_END}"
    echo "  Expanded region: ${EXPANDED_REGION}"

    # Generate PDF with enhanced settings
    echo "  Generating PDF..."
    ggsashimi.py \
        -b ${BAM_CONFIG} \
        -c ${EXPANDED_REGION} \
        -g ${GTF_FILE} \
        -o ${GENE_OUTPUT}/${GENE}_sashimi \
        -M ${MIN_JUNCTION_COV} \
        -C 3 \
        --alpha 0.25 \
        --height ${PLOT_HEIGHT} \
        --width ${PLOT_WIDTH} \
        --base-size 12 \
        --ann-height 2 \
        --shrink \
        -F pdf \
        2>/dev/null

    # Generate PNG with high resolution
    echo "  Generating PNG (${PNG_RESOLUTION} DPI)..."
    ggsashimi.py \
        -b ${BAM_CONFIG} \
        -c ${EXPANDED_REGION} \
        -g ${GTF_FILE} \
        -o ${GENE_OUTPUT}/${GENE}_sashimi \
        -M ${MIN_JUNCTION_COV} \
        -C 3 \
        --alpha 0.25 \
        --height ${PLOT_HEIGHT} \
        --width ${PLOT_WIDTH} \
        --base-size 12 \
        --ann-height 2 \
        --shrink \
        -F png \
        2>/dev/null

    # Verify output
    if [[ -f "${GENE_OUTPUT}/${GENE}_sashimi.pdf" ]]; then
        echo "  SUCCESS: Sashimi plots created"
        ls -lh ${GENE_OUTPUT}/${GENE}_sashimi.*
    else
        echo "  WARNING: Sashimi plot generation may have failed"
    fi
}

#===============================================================================
# Extract Top Differential Splicing Genes
#===============================================================================

extract_top_splicing_genes() {
    echo ""
    echo "=== Extracting top differential splicing genes ==="

    # Python script to extract top genes from rMATS results
    python3 << 'PYEOF'
import os
import pandas as pd
from pathlib import Path

base_dir = os.environ.get('BASE_DIR', '/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649')
splicing_dir = Path(base_dir) / 'results/05_splicing'
output_dir = Path(base_dir) / 'results/06_sashimi'

comparisons = ['EMX1_wt_vs_mut', 'Nestin_wt_vs_mut', 'WT_EMX1_vs_Nestin', 'MUT_EMX1_vs_Nestin']
event_types = ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']

DPSI_THRESHOLD = 0.1
MIN_READS = 20

def sum_counts(count_str):
    if pd.isna(count_str):
        return 0
    try:
        return sum(int(x) for x in str(count_str).split(',') if x.strip())
    except:
        return 0

# Collect top genes
top_genes = []

for comp in comparisons:
    for event in event_types:
        jc_file = splicing_dir / comp / f'{event}.MATS.JC.txt'
        if not jc_file.exists():
            continue

        try:
            df = pd.read_csv(jc_file, sep='\t')

            # Calculate read support
            if all(col in df.columns for col in ['IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2']):
                df['total_reads_1'] = df['IJC_SAMPLE_1'].apply(sum_counts) + df['SJC_SAMPLE_1'].apply(sum_counts)
                df['total_reads_2'] = df['IJC_SAMPLE_2'].apply(sum_counts) + df['SJC_SAMPLE_2'].apply(sum_counts)
                df['min_total_reads'] = df[['total_reads_1', 'total_reads_2']].min(axis=1)
            else:
                df['min_total_reads'] = MIN_READS

            # Filter significant events
            sig = df[
                (abs(df['IncLevelDifference']) > DPSI_THRESHOLD) &
                (df['min_total_reads'] >= MIN_READS)
            ].copy()

            if len(sig) == 0:
                continue

            sig['abs_dpsi'] = abs(sig['IncLevelDifference'])

            # Get top 3 genes by |dPSI|
            top = sig.nlargest(3, 'abs_dpsi')

            for _, row in top.iterrows():
                gene = str(row.get('geneSymbol', row.get('GeneID', 'Unknown'))).replace('"', '')
                if gene and gene != 'nan' and gene != 'Unknown':
                    top_genes.append({
                        'gene': gene,
                        'dpsi': row['IncLevelDifference'],
                        'abs_dpsi': abs(row['IncLevelDifference']),
                        'event': event,
                        'comparison': comp
                    })
        except Exception as e:
            pass

# Sort by |dPSI| and deduplicate (keep highest |dPSI| per gene)
top_genes_df = pd.DataFrame(top_genes)
if len(top_genes_df) > 0:
    top_genes_df = top_genes_df.sort_values('abs_dpsi', ascending=False)
    top_genes_df = top_genes_df.drop_duplicates(subset=['gene'], keep='first')

    # Get top 20 unique genes
    top_20 = top_genes_df.head(20)

    # Save to file
    with open(output_dir / 'top_splicing_genes.txt', 'w') as f:
        for _, row in top_20.iterrows():
            f.write(f"{row['gene']}\n")

    print(f"Extracted {len(top_20)} top differential splicing genes:")
    for _, row in top_20.iterrows():
        print(f"  {row['gene']}: {row['event']} in {row['comparison']}, dPSI={row['dpsi']:.3f}")
else:
    print("No significant splicing events found.")
PYEOF

export BASE_DIR
}

#===============================================================================
# Main Processing
#===============================================================================

# Check if command line arguments provided
if [[ $# -gt 0 ]]; then
    # User-specified genes
    USER_GENES=("$@")
    echo "Processing user-specified genes: ${USER_GENES[@]}"
    echo ""

    for GENE in "${USER_GENES[@]}"; do
        generate_sashimi_enhanced "${GENE}" "" "user_genes"
    done
else
    # Auto-extract and process top splicing genes
    extract_top_splicing_genes

    TOP_GENES_FILE="${OUTPUT_DIR}/top_splicing_genes.txt"
    if [[ -f "${TOP_GENES_FILE}" ]]; then
        echo ""
        echo "Processing top differential splicing genes..."
        while IFS= read -r GENE; do
            if [[ -n "${GENE}" ]]; then
                generate_sashimi_enhanced "${GENE}" "" "top_splicing"
            fi
        done < "${TOP_GENES_FILE}"
    else
        echo "WARNING: No top genes file found. Using default genes."
        # Default housekeeping genes for testing
        DEFAULT_GENES=(
            "Actb"      # Beta-actin
            "Gapdh"     # GAPDH
            "Rbfox1"    # RNA-binding protein
        )
        for GENE in "${DEFAULT_GENES[@]}"; do
            generate_sashimi_enhanced "${GENE}" "" "top_splicing"
        done
    fi
fi

#===============================================================================
# Create Helper Script for Custom Genes
#===============================================================================

HELPER_SCRIPT="${OUTPUT_DIR}/plot_gene_enhanced.sh"
cat > ${HELPER_SCRIPT} << 'HELPER_EOF'
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
HELPER_EOF

chmod +x ${HELPER_SCRIPT}

#===============================================================================
# Summary
#===============================================================================

echo ""
echo "============================================"
echo "Enhanced Sashimi Plot Generation Complete!"
echo "============================================"
echo ""
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "Directories:"
echo "  top_splicing/ - Top differential splicing genes (auto-selected)"
echo "  user_genes/   - User-specified genes"
echo ""
echo "Enhanced settings used:"
echo "  - Plot dimensions: ${PLOT_WIDTH}x${PLOT_HEIGHT} inches"
echo "  - Min junction coverage: ${MIN_JUNCTION_COV}"
echo "  - PNG resolution: ${PNG_RESOLUTION} DPI"
echo "  - Output formats: PDF + PNG"
echo ""
echo "To generate plots for additional genes, use:"
echo "  bash ${HELPER_SCRIPT} GENE_NAME"
echo ""
echo "End time: $(date)"
echo "============================================"
