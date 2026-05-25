#!/bin/bash

#===============================================================================
# SCRIPT: 8_junction_tracks.sh
# PURPOSE: Create IGV-Compatible Junction BED Files from STAR SJ.out.tab
#
# DESCRIPTION:
# Converts STAR splice junction output (SJ.out.tab) files to BED format for
# visualization in IGV and other genome browsers. Filters junctions by minimum
# read count to reduce noise.
#
# INPUT:
# - STAR SJ.out.tab files from alignment step
#   Format: chr, intron_start(1-based), intron_end, strand, motif, annotation_status,
#           unique_reads, multi_reads, overhang
#
# OUTPUT:
# - BED6 files with junction coordinates and read counts as scores
# - Summary statistics for each sample
#
# USAGE:
# sbatch 8_junction_tracks.sh
#===============================================================================

#SBATCH --job-name=8_junction_tracks
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/8_junction_tracks.err"
#SBATCH --output="./logs/8_junction_tracks.out"

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/rnaseq-quant

#===============================================================================
# Configuration
#===============================================================================

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
ALIGNED_DIR="${BASE_DIR}/results/02_aligned"
OUTPUT_DIR="${BASE_DIR}/results/07_igv_junctions"

# Minimum unique reads required to include a junction
MIN_READS=3

# Sample name mapping
declare -A SAMPLE_NAMES
SAMPLE_NAMES[1]="EMX1_hippo_wt"
SAMPLE_NAMES[2]="EMX1_hippo_mut"
SAMPLE_NAMES[3]="Nestin_hippo_wt"
SAMPLE_NAMES[4]="Nestin_hippo_mut"

# Sample IDs
SAMPLES=(1 2 3 4)

#===============================================================================
# Setup
#===============================================================================

echo "==============================================================================="
echo "Creating IGV-Compatible Junction BED Files"
echo "==============================================================================="
echo "Timestamp: $(date)"
echo "Base directory: ${BASE_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Minimum read filter: ${MIN_READS}"
echo ""

# Create output directory
mkdir -p ${OUTPUT_DIR}

#===============================================================================
# Function: Convert SJ.out.tab to BED format
#===============================================================================

convert_sj_to_bed() {
    local SJ_FILE=$1
    local BED_FILE=$2
    local SAMPLE_NAME=$3
    local MIN_READS=$4

    echo "Converting: ${SJ_FILE}"
    echo "Output: ${BED_FILE}"

    # Convert SJ.out.tab to BED6 format
    # SJ.out.tab columns:
    #   1: chromosome
    #   2: intron start (1-based)
    #   3: intron end (1-based, inclusive)
    #   4: strand (0=undefined, 1=+, 2=-)
    #   5: intron motif (0=non-canonical, 1-6=canonical)
    #   6: annotation status (0=novel, 1=annotated)
    #   7: unique mapping reads
    #   8: multi-mapping reads
    #   9: max overhang
    #
    # BED6 output:
    #   1: chromosome
    #   2: start (0-based)
    #   3: end (0-based, exclusive)
    #   4: name (chr:start-end;status)
    #   5: score (unique reads)
    #   6: strand

    awk -v min_reads=${MIN_READS} '
    BEGIN {
        OFS="\t"
        total = 0
        passed = 0
        known = 0
        novel = 0
    }
    {
        total++
        if ($7 >= min_reads) {
            passed++

            # Convert strand: 0=undefined(.), 1=+, 2=-
            if ($4 == 1) {
                strand = "+"
            } else if ($4 == 2) {
                strand = "-"
            } else {
                strand = "."
            }

            # Annotation status: 0=novel, 1=annotated
            if ($6 == 1) {
                annot = "known"
                known++
            } else {
                annot = "novel"
                novel++
            }

            # Create junction name
            name = $1":"$2"-"$3";"annot

            # Convert to 0-based BED coordinates
            # STAR reports 1-based intron start, so subtract 1 for BED
            start = $2 - 1
            end = $3

            # Output BED6 format
            print $1, start, end, name, $7, strand
        }
    }
    END {
        # Print statistics to stderr
        print "  Total junctions: " total > "/dev/stderr"
        print "  Passed filter (>=" min_reads " reads): " passed > "/dev/stderr"
        print "  Known (annotated): " known > "/dev/stderr"
        print "  Novel: " novel > "/dev/stderr"
    }
    ' ${SJ_FILE} | sort -k1,1 -k2,2n > ${BED_FILE}
}

#===============================================================================
# Process each sample
#===============================================================================

echo "==============================================================================="
echo "Processing samples"
echo "==============================================================================="

# Statistics summary file
STATS_FILE="${OUTPUT_DIR}/junction_statistics.txt"
echo -e "Sample\tTotal_Junctions\tFiltered_Junctions\tKnown\tNovel" > ${STATS_FILE}

for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_NAME=${SAMPLE_NAMES[$SAMPLE]}

    echo ""
    echo "--- Processing sample: ${SAMPLE} (${SAMPLE_NAME}) ---"

    # Input SJ.out.tab file
    SJ_FILE="${ALIGNED_DIR}/${SAMPLE}/${SAMPLE}_SJ.out.tab"

    # Output BED file
    BED_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_junctions.bed"

    # Check if input file exists
    if [[ ! -f "${SJ_FILE}" ]]; then
        echo "WARNING: SJ.out.tab file not found: ${SJ_FILE}"
        echo "  Skipping sample ${SAMPLE}"
        continue
    fi

    # Convert to BED format
    convert_sj_to_bed ${SJ_FILE} ${BED_FILE} ${SAMPLE_NAME} ${MIN_READS}

    # Verify output
    if [[ -f "${BED_FILE}" ]]; then
        JUNCTION_COUNT=$(wc -l < ${BED_FILE})
        echo "  Output file: ${BED_FILE}"
        echo "  Junction count: ${JUNCTION_COUNT}"

        # Add to statistics
        TOTAL=$(wc -l < ${SJ_FILE})
        KNOWN=$(grep -c ";known$" ${BED_FILE} || echo 0)
        NOVEL=$(grep -c ";novel$" ${BED_FILE} || echo 0)
        echo -e "${SAMPLE_NAME}\t${TOTAL}\t${JUNCTION_COUNT}\t${KNOWN}\t${NOVEL}" >> ${STATS_FILE}
    else
        echo "  ERROR: Failed to create BED file"
    fi
done

#===============================================================================
# Create merged junction file from all samples
#===============================================================================

echo ""
echo "==============================================================================="
echo "Creating merged junction file"
echo "==============================================================================="

MERGED_FILE="${OUTPUT_DIR}/all_samples_merged_junctions.bed"

# Collect all junction coordinates, merge overlapping junctions, and sum read counts
# We'll use a more sophisticated approach: track junction positions and combine scores

echo "Merging junctions from all samples..."

# Create a temporary file with all junctions labeled by sample
TEMP_ALL="${OUTPUT_DIR}/temp_all_junctions.txt"
rm -f ${TEMP_ALL}

for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_NAME=${SAMPLE_NAMES[$SAMPLE]}
    BED_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_junctions.bed"

    if [[ -f "${BED_FILE}" ]]; then
        # Add sample name as column 7
        awk -v sample=${SAMPLE_NAME} 'BEGIN{OFS="\t"} {print $0, sample}' ${BED_FILE} >> ${TEMP_ALL}
    fi
done

# Merge junctions: combine by exact coordinates, sum scores, list samples
echo "Aggregating junction data..."

awk '
BEGIN {
    OFS="\t"
}
{
    # Key by chr, start, end, strand
    key = $1"\t"$2"\t"$3"\t"$6

    # Sum the scores (read counts)
    scores[key] += $5

    # Track samples that have this junction
    if (samples[key] == "") {
        samples[key] = $7
    } else if (index(samples[key], $7) == 0) {
        samples[key] = samples[key]","$7
    }

    # Keep track of annotation status (prefer "known" if any sample has it)
    # Extract status from name field ($4 which is chr:start-end;status)
    split($4, parts, ";")
    if (parts[2] == "known") {
        known[key] = 1
    }

    # Store original name for reconstruction
    chr[key] = $1
    start[key] = $2
    end_pos[key] = $3
    strand[key] = $6
}
END {
    for (key in scores) {
        status = (known[key] == 1) ? "known" : "novel"
        name = chr[key]":"(start[key]+1)"-"end_pos[key]";"status";n="length(split(samples[key], a, ","))
        print chr[key], start[key], end_pos[key], name, scores[key], strand[key]
    }
}
' ${TEMP_ALL} | sort -k1,1 -k2,2n > ${MERGED_FILE}

# Cleanup
rm -f ${TEMP_ALL}

# Report merged stats
MERGED_COUNT=$(wc -l < ${MERGED_FILE})
echo "  Merged junction file: ${MERGED_FILE}"
echo "  Total unique junctions: ${MERGED_COUNT}"

#===============================================================================
# Summary
#===============================================================================

echo ""
echo "==============================================================================="
echo "Summary"
echo "==============================================================================="

echo ""
echo "Junction statistics:"
cat ${STATS_FILE}

echo ""
echo "Generated files:"
ls -lh ${OUTPUT_DIR}/*.bed 2>/dev/null || echo "  No BED files generated"

echo ""
echo "==============================================================================="
echo "Junction track generation complete"
echo "==============================================================================="
echo "Timestamp: $(date)"
echo ""
echo "To load in IGV:"
echo "  1. Open IGV and load the mouse genome (mm39/GRCm39)"
echo "  2. File > Load from File > Select the .bed files"
echo "  3. Junction tracks will display with read counts as scores"
echo "  4. Use 'Sashimi Plot' view for detailed junction visualization"
echo ""
echo "Files are located in: ${OUTPUT_DIR}"
