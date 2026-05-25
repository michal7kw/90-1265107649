#!/bin/bash
#SBATCH --job-name=12_igv_subset
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/12_igv_subset.err"
#SBATCH --output="./logs/12_igv_subset.out"

#===============================================================================
# Script: 12_create_igv_subset.sh
# Purpose: Create a lightweight IGV package with subsetted BAM files
#
# Description:
# This script extracts only the genomic regions containing significant
# differential splicing events from the full BAM files. This dramatically
# reduces file size while preserving the ability to view sashimi plots
# for all significant events.
#
# Output:
# - Subsetted BAM files (~100-500 MB instead of ~8 GB each)
# - BigWig coverage tracks
# - BED files with splicing events
# - IGV session XML files
# - Coordinate navigation files
#
# Usage:
# bash scripts/12_create_igv_subset.sh
# or
# sbatch scripts/12_create_igv_subset.sh
#===============================================================================

set -euo pipefail

echo "============================================"
echo "Creating IGV Subset Package"
echo "Start time: $(date)"
echo "============================================"

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
ALIGNED_DIR="${BASE_DIR}/results/02_aligned"
IGV_DIR="${BASE_DIR}/results/08_igv"
BIGWIG_DIR="${BASE_DIR}/results/03_bigwig"
OUTPUT_DIR="${BASE_DIR}/results/08_igv_subset"

# Sample mapping
declare -A SAMPLES
SAMPLES[1]="EMX1_wt"
SAMPLES[2]="EMX1_mut"
SAMPLES[3]="Nestin_wt"
SAMPLES[4]="Nestin_mut"

# Padding around splicing events (bp)
PADDING=2000

#-------------------------------------------------------------------------------
# Setup conda environment
#-------------------------------------------------------------------------------

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/ggsashimi

# Verify samtools is available
if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found in PATH"
    exit 1
fi

echo "Using samtools: $(which samtools)"
echo "Samtools version: $(samtools --version | head -1)"
echo ""

#-------------------------------------------------------------------------------
# Create output directory structure
#-------------------------------------------------------------------------------

echo "=== Creating output directory structure ==="
mkdir -p ${OUTPUT_DIR}/{bam,bed,bigwig,coordinates}

#-------------------------------------------------------------------------------
# Step 1: Create merged regions BED file from all significant splicing events
#-------------------------------------------------------------------------------

echo ""
echo "=== Step 1: Creating merged regions BED file ==="

# Extract coordinates from all significant splicing BED files
# Add padding and merge overlapping regions
cat ${IGV_DIR}/bed/*_significant_splicing.bed | \
    grep -v "^track" | \
    awk -v padding=${PADDING} -v OFS='\t' '{
        # Add padding on each side for visualization context
        start = ($2 - padding < 0) ? 0 : $2 - padding;
        end = $3 + padding;
        print $1, start, end
    }' | \
    sort -k1,1 -k2,2n | \
    uniq > ${OUTPUT_DIR}/splicing_regions_unsorted.bed

# Merge overlapping regions using awk (since bedtools may not be available)
# This is a simple merge for sorted BED files
awk -v OFS='\t' '
    NR == 1 {
        chr = $1; start = $2; end = $3; next
    }
    $1 == chr && $2 <= end {
        # Overlapping or adjacent - extend
        if ($3 > end) end = $3
    }
    $1 != chr || $2 > end {
        # New region - print previous and start new
        print chr, start, end
        chr = $1; start = $2; end = $3
    }
    END {
        print chr, start, end
    }
' ${OUTPUT_DIR}/splicing_regions_unsorted.bed > ${OUTPUT_DIR}/splicing_regions.bed

# Cleanup
rm -f ${OUTPUT_DIR}/splicing_regions_unsorted.bed

# Report statistics
N_REGIONS=$(wc -l < ${OUTPUT_DIR}/splicing_regions.bed)
TOTAL_BP=$(awk '{sum += $3-$2} END {print sum}' ${OUTPUT_DIR}/splicing_regions.bed)
TOTAL_MB=$(echo "scale=2; ${TOTAL_BP}/1000000" | bc)

echo "  Merged regions: ${N_REGIONS}"
echo "  Total bases: ${TOTAL_MB} Mb"
echo "  Regions file: ${OUTPUT_DIR}/splicing_regions.bed"

#-------------------------------------------------------------------------------
# Step 2: Extract BAM subsets for each sample
#-------------------------------------------------------------------------------

echo ""
echo "=== Step 2: Extracting BAM subsets ==="

for SAMPLE_ID in "${!SAMPLES[@]}"; do
    SAMPLE_NAME="${SAMPLES[$SAMPLE_ID]}"
    INPUT_BAM="${ALIGNED_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam"
    OUTPUT_BAM="${OUTPUT_DIR}/bam/${SAMPLE_NAME}.bam"

    echo ""
    echo "Processing ${SAMPLE_NAME} (Sample ${SAMPLE_ID})..."

    if [[ ! -f "${INPUT_BAM}" ]]; then
        echo "  WARNING: Input BAM not found: ${INPUT_BAM}"
        continue
    fi

    # Extract reads from regions of interest
    echo "  Extracting reads from ${N_REGIONS} regions..."
    samtools view -b -h -L ${OUTPUT_DIR}/splicing_regions.bed \
        -@ 4 \
        ${INPUT_BAM} > ${OUTPUT_BAM}

    # Index the subsetted BAM
    echo "  Indexing..."
    samtools index -@ 4 ${OUTPUT_BAM}

    # Report sizes
    ORIG_SIZE=$(du -h ${INPUT_BAM} | cut -f1)
    NEW_SIZE=$(du -h ${OUTPUT_BAM} | cut -f1)
    echo "  Original: ${ORIG_SIZE} → Subset: ${NEW_SIZE}"
done

#-------------------------------------------------------------------------------
# Step 3: Copy BigWig files
#-------------------------------------------------------------------------------

echo ""
echo "=== Step 3: Copying BigWig files ==="

for SAMPLE_ID in "${!SAMPLES[@]}"; do
    SAMPLE_NAME="${SAMPLES[$SAMPLE_ID]}"

    # Try different naming conventions
    if [[ -f "${BIGWIG_DIR}/${SAMPLE_NAME}_RPKM.bw" ]]; then
        SRC="${BIGWIG_DIR}/${SAMPLE_NAME}_RPKM.bw"
    elif [[ -f "${BIGWIG_DIR}/${SAMPLE_NAME//_/_hippo_}_RPKM.bw" ]]; then
        # Handle EMX1_wt -> EMX1_hippo_wt naming
        HIPPO_NAME=$(echo ${SAMPLE_NAME} | sed 's/_/_hippo_/')
        SRC="${BIGWIG_DIR}/${HIPPO_NAME}_RPKM.bw"
    else
        echo "  WARNING: BigWig not found for ${SAMPLE_NAME}"
        continue
    fi

    DST="${OUTPUT_DIR}/bigwig/${SAMPLE_NAME}.bw"
    cp "${SRC}" "${DST}"
    echo "  Copied: ${SAMPLE_NAME}.bw ($(du -h ${DST} | cut -f1))"
done

#-------------------------------------------------------------------------------
# Step 4: Copy BED files and coordinates
#-------------------------------------------------------------------------------

echo ""
echo "=== Step 4: Copying BED files and coordinates ==="

# Copy BED files
cp ${IGV_DIR}/bed/*.bed ${OUTPUT_DIR}/bed/
echo "  Copied $(ls ${OUTPUT_DIR}/bed/*.bed | wc -l) BED files"

# Copy coordinate files
cp ${IGV_DIR}/coordinates/* ${OUTPUT_DIR}/coordinates/ 2>/dev/null || true
cp ${IGV_DIR}/top_splicing_regions.txt ${OUTPUT_DIR}/ 2>/dev/null || true
echo "  Copied coordinate files"

#-------------------------------------------------------------------------------
# Step 5: Create IGV session XML with relative paths
#-------------------------------------------------------------------------------

echo ""
echo "=== Step 5: Creating IGV session files ==="

# Full session with BAM files (for sashimi plots)
cat > ${OUTPUT_DIR}/igv_session_with_bam.xml << 'XMLEOF'
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="mm39" hasGeneTrack="true" hasSequenceTrack="true" locus="All" version="8">
    <Resources>
        <!-- BigWig Coverage Tracks -->
        <Resource path="bigwig/EMX1_wt.bw"/>
        <Resource path="bigwig/EMX1_mut.bw"/>
        <Resource path="bigwig/Nestin_wt.bw"/>
        <Resource path="bigwig/Nestin_mut.bw"/>
        <!-- BAM files (for Sashimi plots) -->
        <Resource path="bam/EMX1_wt.bam"/>
        <Resource path="bam/EMX1_mut.bam"/>
        <Resource path="bam/Nestin_wt.bam"/>
        <Resource path="bam/Nestin_mut.bam"/>
        <!-- Splicing Event BED Tracks -->
        <Resource path="bed/EMX1_wt_vs_mut_significant_splicing.bed"/>
        <Resource path="bed/Nestin_wt_vs_mut_significant_splicing.bed"/>
        <Resource path="bed/WT_EMX1_vs_Nestin_significant_splicing.bed"/>
        <Resource path="bed/MUT_EMX1_vs_Nestin_significant_splicing.bed"/>
    </Resources>
    <Panel name="DataPanel" height="400">
    </Panel>
    <Panel name="FeaturePanel" height="150">
    </Panel>
    <PanelLayout dividerFractions="0.72"/>
</Session>
XMLEOF

# Light session without BAM files (just coverage and events)
cat > ${OUTPUT_DIR}/igv_session_light.xml << 'XMLEOF'
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="mm39" hasGeneTrack="true" hasSequenceTrack="true" locus="All" version="8">
    <Resources>
        <!-- BigWig Coverage Tracks -->
        <Resource path="bigwig/EMX1_wt.bw"/>
        <Resource path="bigwig/EMX1_mut.bw"/>
        <Resource path="bigwig/Nestin_wt.bw"/>
        <Resource path="bigwig/Nestin_mut.bw"/>
        <!-- Splicing Event BED Tracks -->
        <Resource path="bed/EMX1_wt_vs_mut_significant_splicing.bed"/>
        <Resource path="bed/Nestin_wt_vs_mut_significant_splicing.bed"/>
        <Resource path="bed/WT_EMX1_vs_Nestin_significant_splicing.bed"/>
        <Resource path="bed/MUT_EMX1_vs_Nestin_significant_splicing.bed"/>
    </Resources>
    <Panel name="DataPanel" height="400">
    </Panel>
    <Panel name="FeaturePanel" height="150">
    </Panel>
    <PanelLayout dividerFractions="0.72"/>
</Session>
XMLEOF

echo "  Created: igv_session_with_bam.xml"
echo "  Created: igv_session_light.xml"

#-------------------------------------------------------------------------------
# Step 6: Create README file
#-------------------------------------------------------------------------------

echo ""
echo "=== Step 6: Creating README ==="

cat > ${OUTPUT_DIR}/README.txt << EOF
================================================================================
IGV Visualization Subset Package
Project: 90-1265107649
Created: $(date)
================================================================================

DESCRIPTION
-----------
This package contains subsetted BAM files that include ONLY the genomic regions
with significant differential splicing events. This reduces file size from
~35 GB to ~1-2 GB while preserving full sashimi plot capability for all
significant events.

CONTENTS
--------
bam/                    - Subsetted BAM files + indexes (for sashimi plots)
bigwig/                 - Coverage tracks (RPKM normalized)
bed/                    - Color-coded splicing event BED files
coordinates/            - IGV navigation coordinate files
splicing_regions.bed    - Regions extracted from BAMs
igv_session_with_bam.xml - IGV session with BAM files (for sashimi)
igv_session_light.xml   - IGV session without BAMs (lighter)
top_splicing_regions.txt - Quick navigation file

SAMPLE MAPPING
--------------
EMX1_wt.bam    = Sample 1 - EMX1 promoter, wild-type
EMX1_mut.bam   = Sample 2 - EMX1 promoter, mutant
Nestin_wt.bam  = Sample 3 - Nestin promoter, wild-type
Nestin_mut.bam = Sample 4 - Nestin promoter, mutant

COLOR CODING (BED TRACKS)
-------------------------
SE   = Red (255,0,0)       - Skipped Exon
A5SS = Green (0,255,0)     - Alternative 5' Splice Site
A3SS = Blue (0,0,255)      - Alternative 3' Splice Site
MXE  = Orange (255,165,0)  - Mutually Exclusive Exons
RI   = Purple (128,0,128)  - Retained Intron

USAGE
-----
1. Download this entire folder to your local machine
2. Open IGV and set genome to mm39 (GRCm39)
3. File > Load Session > igv_session_with_bam.xml
4. Navigate using coordinates from coordinates/IGV_coordinates_all.txt
5. Right-click on BAM track > Sashimi Plot

STATISTICS
----------
Regions extracted: ${N_REGIONS}
Total bases covered: ${TOTAL_MB} Mb
Padding per region: ${PADDING} bp

IMPORTANT NOTE
--------------
The subsetted BAM files only contain reads from significant splicing regions.
If you navigate to a region NOT in splicing_regions.bed, the BAM tracks will
appear empty. Use the coordinate files to navigate to valid regions.

================================================================================
EOF

echo "  Created: README.txt"

#-------------------------------------------------------------------------------
# Step 7: Create tarball for easy download
#-------------------------------------------------------------------------------

echo ""
echo "=== Step 7: Creating download tarball ==="

cd ${BASE_DIR}/results
tar -cvzf 08_igv_subset.tar.gz 08_igv_subset/

TARBALL_SIZE=$(du -h 08_igv_subset.tar.gz | cut -f1)
echo ""
echo "  Tarball created: 08_igv_subset.tar.gz (${TARBALL_SIZE})"

#-------------------------------------------------------------------------------
# Summary
#-------------------------------------------------------------------------------

echo ""
echo "============================================"
echo "IGV Subset Package Complete!"
echo "============================================"
echo ""
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "Package contents:"
du -sh ${OUTPUT_DIR}/*
echo ""
echo "Total package size:"
du -sh ${OUTPUT_DIR}
echo ""
echo "Tarball for download:"
echo "  ${BASE_DIR}/results/08_igv_subset.tar.gz (${TARBALL_SIZE})"
echo ""
echo "Download command:"
echo "  rsync -avP kubacki.michal@srhpclogin01.ihsr.dom:${BASE_DIR}/results/08_igv_subset.tar.gz ./"
echo ""
echo "End time: $(date)"
echo "============================================"
