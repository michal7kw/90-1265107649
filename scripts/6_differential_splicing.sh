#!/bin/bash

#===============================================================================
# SCRIPT: 6_differential_splicing.sh
# PURPOSE: Differential Splicing Analysis using rMATS
#
# DESCRIPTION:
# Performs differential alternative splicing analysis between conditions
# using rMATS (replicate Multivariate Analysis of Transcript Splicing).
#
# SPLICING EVENTS DETECTED:
# - SE: Skipped Exon
# - A5SS: Alternative 5' Splice Site
# - A3SS: Alternative 3' Splice Site
# - MXE: Mutually Exclusive Exons
# - RI: Retained Intron
#
# COMPARISONS:
# 1. EMX1: wt vs mut
# 2. Nestin: wt vs mut
# 3. WT: EMX1 vs Nestin
# 4. MUT: EMX1 vs Nestin
#
# IMPORTANT - MONOPLICATE DESIGN LIMITATIONS:
# This experiment has n=1 per condition (no biological replicates).
# rMATS is designed for replicated experiments and uses statistical tests.
# With monoplicates:
#   - P-values are NOT reliable (set tstat=1 to minimize but not eliminate)
#   - High number of "significant" events expected (mostly noise)
#   - Focus ONLY on PSI differences > 0.1 (10%) with high read support
#   - Results should be treated as EXPLORATORY, not definitive
#   - Validation of candidate events is REQUIRED
#
# LIBRARY STRANDEDNESS:
# Verified as UNSTRANDED (Forward/Reverse ratio ≈ 0.97)
# Using --libType fr-unstranded
#
# USAGE:
# sbatch 6_differential_splicing.sh
#===============================================================================

#SBATCH --job-name=6_splicing_ribotag
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/6_splicing.err"
#SBATCH --output="./logs/6_splicing.out"

# Load conda
source /opt/common/tools/ric.cosr/miniconda3/bin/activate

# Check if rmats environment exists, if not create it
RMATS_ENV="/beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/rmats"

if [[ ! -d "${RMATS_ENV}" ]]; then
    echo "Creating rMATS conda environment..."
    conda create -y -p ${RMATS_ENV} -c bioconda -c conda-forge rmats=4.3.0 python=3.9
fi

conda activate ${RMATS_ENV}

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
ALIGNED_DIR="${BASE_DIR}/results/02_aligned"
OUTPUT_DIR="${BASE_DIR}/results/05_splicing"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "=== Starting Differential Splicing Analysis ==="
echo "Timestamp: $(date)"

#-------------------------------------------------------------------------------
# Function to run rMATS comparison
#-------------------------------------------------------------------------------
run_rmats_comparison() {
    local SAMPLE1=$1
    local SAMPLE2=$2
    local COMPARISON_NAME=$3

    echo ""
    echo "=== Running comparison: ${COMPARISON_NAME} ==="
    echo "Sample 1: ${SAMPLE1}"
    echo "Sample 2: ${SAMPLE2}"

    # Create BAM list files for rMATS
    local COMP_DIR="${OUTPUT_DIR}/${COMPARISON_NAME}"
    mkdir -p ${COMP_DIR}

    # Get BAM file paths
    local BAM1="${ALIGNED_DIR}/${SAMPLE1}/${SAMPLE1}_Aligned.sortedByCoord.out.bam"
    local BAM2="${ALIGNED_DIR}/${SAMPLE2}/${SAMPLE2}_Aligned.sortedByCoord.out.bam"

    # Create b1.txt and b2.txt files
    echo "${BAM1}" > ${COMP_DIR}/b1.txt
    echo "${BAM2}" > ${COMP_DIR}/b2.txt

    # Run rMATS
    # --readLength should match your read length (150bp based on earlier check)
    # --variable-read-length for variable length reads
    # -t paired for paired-end data
    rmats.py \
        --b1 ${COMP_DIR}/b1.txt \
        --b2 ${COMP_DIR}/b2.txt \
        --gtf ${GTF_FILE} \
        --od ${COMP_DIR} \
        --tmp ${COMP_DIR}/tmp \
        -t paired \
        --readLength 150 \
        --variable-read-length \
        --nthread 16 \
        --tstat 1 \
        --cstat 0.1 \
        --libType fr-unstranded

    # Check if completed
    if [[ -f "${COMP_DIR}/SE.MATS.JC.txt" ]]; then
        echo "SUCCESS: rMATS completed for ${COMPARISON_NAME}"

        # Count significant events with stricter filtering for monoplicates
        # Require: |deltaPSI| > 0.1 AND minimum 20 junction reads
        for EVENT in SE A5SS A3SS MXE RI; do
            if [[ -f "${COMP_DIR}/${EVENT}.MATS.JC.txt" ]]; then
                # Count events with |deltaPSI| > 0.1 AND min reads >= 20
                # Columns: 13=IJC_S1, 14=SJC_S1, 15=IJC_S2, 16=SJC_S2, 23=IncLevelDiff
                COUNT=$(awk -F'\t' 'NR>1 && ($23 > 0.1 || $23 < -0.1) && ($13+$14 >= 20) && ($15+$16 >= 20) {count++} END {print count+0}' \
                    ${COMP_DIR}/${EVENT}.MATS.JC.txt)
                echo "  ${EVENT}: ${COUNT} events (|deltaPSI| > 0.1, min reads >= 20)"
            fi
        done
    else
        echo "WARNING: rMATS may have failed for ${COMPARISON_NAME}"
    fi

    # Clean up tmp directory
    rm -rf ${COMP_DIR}/tmp
}

#-------------------------------------------------------------------------------
# Run all comparisons
#-------------------------------------------------------------------------------

# Comparison 1: EMX1 wt vs mut
run_rmats_comparison "1" "2" "EMX1_wt_vs_mut"

# Comparison 2: Nestin wt vs mut
run_rmats_comparison "3" "4" "Nestin_wt_vs_mut"

# Comparison 3: WT - EMX1 vs Nestin
run_rmats_comparison "1" "3" "WT_EMX1_vs_Nestin"

# Comparison 4: MUT - EMX1 vs Nestin
run_rmats_comparison "2" "4" "MUT_EMX1_vs_Nestin"

#-------------------------------------------------------------------------------
# Generate summary
#-------------------------------------------------------------------------------
echo ""
echo "=== Generating Summary ==="

SUMMARY_FILE="${OUTPUT_DIR}/splicing_summary.txt"
echo "Differential Splicing Analysis Summary" > ${SUMMARY_FILE}
echo "======================================" >> ${SUMMARY_FILE}
echo "Date: $(date)" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}
echo "IMPORTANT: Monoplicate data (n=1) - results are EXPLORATORY only" >> ${SUMMARY_FILE}
echo "Filtering: |deltaPSI| > 0.1 AND min junction reads >= 20" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}

for COMP in EMX1_wt_vs_mut Nestin_wt_vs_mut WT_EMX1_vs_Nestin MUT_EMX1_vs_Nestin; do
    echo "Comparison: ${COMP}" >> ${SUMMARY_FILE}
    echo "------------------------" >> ${SUMMARY_FILE}
    for EVENT in SE A5SS A3SS MXE RI; do
        if [[ -f "${OUTPUT_DIR}/${COMP}/${EVENT}.MATS.JC.txt" ]]; then
            # Stricter filtering for monoplicate data
            COUNT=$(awk -F'\t' 'NR>1 && ($23 > 0.1 || $23 < -0.1) && ($13+$14 >= 20) && ($15+$16 >= 20) {count++} END {print count+0}' \
                ${OUTPUT_DIR}/${COMP}/${EVENT}.MATS.JC.txt)
            echo "  ${EVENT}: ${COUNT}" >> ${SUMMARY_FILE}
        fi
    done
    echo "" >> ${SUMMARY_FILE}
done

cat ${SUMMARY_FILE}

echo ""
echo "=== Differential Splicing Analysis Complete ==="
echo "Results in: ${OUTPUT_DIR}/"
echo "Timestamp: $(date)"
