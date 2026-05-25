#!/bin/bash

#===============================================================================
# SCRIPT: check_strandedness.sh
# PURPOSE: Detect library strandedness from STAR gene counts
#
# DESCRIPTION:
# Analyzes STAR ReadsPerGene.out.tab files to determine library strandedness.
# - Unstranded: Forward/Reverse ratio ≈ 1.0 (0.8-1.2)
# - fr-firststrand (dUTP/TruSeq): Forward/Reverse ratio < 0.3
# - fr-secondstrand: Forward/Reverse ratio > 3.0
#
# USAGE:
# bash check_strandedness.sh
#===============================================================================

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
ALIGNED_DIR="${BASE_DIR}/results/02_aligned"

echo "=============================================="
echo "Library Strandedness Analysis"
echo "=============================================="
echo ""
echo "Interpretation guide:"
echo "  - Ratio ≈ 1.0 (0.8-1.2): UNSTRANDED"
echo "  - Ratio < 0.3: fr-firststrand (e.g., dUTP, TruSeq stranded)"
echo "  - Ratio > 3.0: fr-secondstrand"
echo ""
echo "----------------------------------------------"

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))

# Arrays to store results
declare -a RATIOS

for SAMPLE in "${SAMPLES[@]}"; do
    COUNTS_FILE="${ALIGNED_DIR}/${SAMPLE}/${SAMPLE}_ReadsPerGene.out.tab"

    if [[ ! -f "${COUNTS_FILE}" ]]; then
        echo "WARNING: Counts file not found for sample ${SAMPLE}"
        continue
    fi

    # Calculate sums and ratio
    RESULT=$(awk 'NR>4 {sum2+=$2; sum3+=$3; sum4+=$4} END {
        ratio = sum3/sum4
        printf "%d\t%d\t%d\t%.4f", sum2, sum3, sum4, ratio
    }' ${COUNTS_FILE})

    UNSTRANDED=$(echo "$RESULT" | cut -f1)
    FORWARD=$(echo "$RESULT" | cut -f2)
    REVERSE=$(echo "$RESULT" | cut -f3)
    RATIO=$(echo "$RESULT" | cut -f4)

    RATIOS+=($RATIO)

    printf "Sample %s:\n" "${SAMPLE}"
    printf "  Unstranded counts: %'d\n" "${UNSTRANDED}"
    printf "  Forward counts:    %'d\n" "${FORWARD}"
    printf "  Reverse counts:    %'d\n" "${REVERSE}"
    printf "  Forward/Reverse:   %.4f\n" "${RATIO}"
    echo ""
done

echo "----------------------------------------------"
echo "SUMMARY"
echo "----------------------------------------------"

# Calculate average ratio
AVG_RATIO=$(echo "${RATIOS[@]}" | tr ' ' '\n' | awk '{sum+=$1; count++} END {print sum/count}')

printf "Average Forward/Reverse ratio: %.4f\n" "${AVG_RATIO}"
echo ""

# Determine strandedness
if (( $(echo "$AVG_RATIO > 0.8 && $AVG_RATIO < 1.2" | bc -l) )); then
    echo "CONCLUSION: Library is UNSTRANDED"
    echo ""
    echo "Recommended settings:"
    echo "  - STAR counts: Use column 2 (unstranded)"
    echo "  - rMATS: --libType fr-unstranded"
    echo "  - featureCounts: -s 0"
elif (( $(echo "$AVG_RATIO < 0.3" | bc -l) )); then
    echo "CONCLUSION: Library is fr-firststrand (e.g., dUTP/TruSeq stranded)"
    echo ""
    echo "Recommended settings:"
    echo "  - STAR counts: Use column 4 (reverse)"
    echo "  - rMATS: --libType fr-firststrand"
    echo "  - featureCounts: -s 2"
elif (( $(echo "$AVG_RATIO > 3.0" | bc -l) )); then
    echo "CONCLUSION: Library is fr-secondstrand"
    echo ""
    echo "Recommended settings:"
    echo "  - STAR counts: Use column 3 (forward)"
    echo "  - rMATS: --libType fr-secondstrand"
    echo "  - featureCounts: -s 1"
else
    echo "CONCLUSION: Strandedness UNCLEAR (ratio between 0.3-0.8 or 1.2-3.0)"
    echo "Consider using RSeQC infer_experiment.py for more accurate detection"
fi

echo ""
echo "=============================================="
