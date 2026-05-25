#!/bin/bash
#SBATCH --job-name=10_igv_prep
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/10_igv_prep.err"
#SBATCH --output="./logs/10_igv_prep.out"

# =============================================================================
# Script: 10_prepare_igv.sh
# Purpose: Prepare files for IGV visualization of differential splicing
#
# Creates:
#   1. Color-coded BED files from rMATS results (by event type)
#   2. IGV session XML file
#   3. Top splicing regions navigation file
#
# IMPORTANT - MONOPLICATE DESIGN:
# This experiment has n=1 per condition. Filtering uses:
#   - |deltaPSI| > 0.1 (10% difference threshold)
#   - Minimum junction reads >= 20
# Results should be treated as exploratory.
# =============================================================================

set -euo pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
BAM_DIR="${BASE_DIR}/results/02_aligned"
BIGWIG_DIR="${BASE_DIR}/results/03_bigwig"
SPLICING_DIR="${BASE_DIR}/results/05_splicing"
OUTPUT_DIR="${BASE_DIR}/results/08_igv"

echo "============================================"
echo "Preparing IGV Visualization Files"
echo "Start time: $(date)"
echo "============================================"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/rmats

# Create output directories
mkdir -p ${OUTPUT_DIR}/{bed,bigwig}

# =============================================================================
# 1. Create Color-Coded BED Files from rMATS Results
# =============================================================================
echo ""
echo "=== Creating BED files for splicing events ==="

python3 << 'PYEOF'
import os
import pandas as pd
from pathlib import Path

base_dir = os.environ.get('BASE_DIR', '/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649')
splicing_dir = Path(base_dir) / 'results/05_splicing'
output_dir = Path(base_dir) / 'results/08_igv/bed'

# Comparisons for this project (monoplicate design)
comparisons = ['EMX1_wt_vs_mut', 'Nestin_wt_vs_mut', 'WT_EMX1_vs_Nestin', 'MUT_EMX1_vs_Nestin']
event_types = ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']

# Color scheme for different event types (RGB)
# These colors will be visible in IGV when itemRgb is enabled
colors = {
    'SE': '255,0,0',      # Red - Skipped Exon
    'A5SS': '0,255,0',    # Green - Alt 5' splice site
    'A3SS': '0,0,255',    # Blue - Alt 3' splice site
    'MXE': '255,165,0',   # Orange - Mutually exclusive exons
    'RI': '128,0,128'     # Purple - Retained intron
}

# Filtering thresholds for monoplicate data (no FDR available)
DPSI_THRESHOLD = 0.1
MIN_READS = 20

def sum_counts(count_str):
    """Sum comma-separated count values."""
    if pd.isna(count_str):
        return 0
    try:
        return sum(int(x) for x in str(count_str).split(',') if x.strip())
    except:
        return 0

for comp in comparisons:
    all_events = []

    for event in event_types:
        jc_file = splicing_dir / comp / f'{event}.MATS.JC.txt'
        if not jc_file.exists():
            print(f"  Skipping {comp}/{event} (file not found)")
            continue

        try:
            df = pd.read_csv(jc_file, sep='\t')

            # Calculate total junction reads for filtering
            if all(col in df.columns for col in ['IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2']):
                df['total_reads_1'] = df['IJC_SAMPLE_1'].apply(sum_counts) + df['SJC_SAMPLE_1'].apply(sum_counts)
                df['total_reads_2'] = df['IJC_SAMPLE_2'].apply(sum_counts) + df['SJC_SAMPLE_2'].apply(sum_counts)
                df['min_total_reads'] = df[['total_reads_1', 'total_reads_2']].min(axis=1)
            else:
                df['min_total_reads'] = MIN_READS  # Default to threshold if columns missing

            # Filter significant events using deltaPSI threshold (not FDR for monoplicate)
            # Also require minimum read support
            sig = df[
                (abs(df['IncLevelDifference']) > DPSI_THRESHOLD) &
                (df['min_total_reads'] >= MIN_READS)
            ].copy()

            if len(sig) == 0:
                continue

            for _, row in sig.iterrows():
                chrom = row['chr']

                # Get coordinates based on event type
                if event == 'SE':
                    start = int(row['exonStart_0base'])
                    end = int(row['exonEnd'])
                elif event == 'RI':
                    start = int(row['riExonStart_0base'])
                    end = int(row['riExonEnd'])
                elif event in ['A5SS', 'A3SS']:
                    start = min(int(row['longExonStart_0base']), int(row['shortES']))
                    end = max(int(row['longExonEnd']), int(row['shortEE']))
                elif event == 'MXE':
                    start = min(int(row['1stExonStart_0base']), int(row['2ndExonStart_0base']))
                    end = max(int(row['1stExonEnd']), int(row['2ndExonEnd']))

                gene = row.get('geneSymbol', row.get('GeneID', 'Unknown'))
                if isinstance(gene, str):
                    gene = gene.replace('"', '')

                dpsi = row['IncLevelDifference']

                # BED9 format: chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb
                name = f"{gene}_{event}_dPSI={dpsi:.2f}"
                score = min(1000, int(abs(dpsi) * 1000))  # Scale |dPSI| to score (0-1000)
                strand = row.get('strand', '.')

                all_events.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': name,
                    'score': score,
                    'strand': strand,
                    'thickStart': start,
                    'thickEnd': end,
                    'itemRgb': colors[event],
                    'event_type': event,
                    'dpsi': dpsi,
                    'gene': gene
                })
        except Exception as e:
            print(f"  Error processing {jc_file}: {e}")

    if all_events:
        bed_df = pd.DataFrame(all_events)

        # Write combined BED file for this comparison (BED9 format with colors)
        bed_file = output_dir / f'{comp}_significant_splicing.bed'
        with open(bed_file, 'w') as f:
            f.write(f'track name="{comp}_splicing" description="Splicing events (|dPSI|>{DPSI_THRESHOLD}, reads>={MIN_READS})" itemRgb="On"\n')
        bed_df[['chrom', 'start', 'end', 'name', 'score', 'strand',
                'thickStart', 'thickEnd', 'itemRgb']].to_csv(
            bed_file, sep='\t', header=False, index=False, mode='a'
        )
        print(f"  {comp}: {len(all_events)} significant events written to BED")

        # Also write separate BED files per event type
        for event in event_types:
            event_df = bed_df[bed_df['event_type'] == event]
            if len(event_df) > 0:
                event_file = output_dir / f'{comp}_{event}.bed'
                rgb_space = colors[event].replace(',', ' ')
                with open(event_file, 'w') as f:
                    f.write(f'track name="{comp}_{event}" description="{event} events" color="{rgb_space}"\n')
                event_df[['chrom', 'start', 'end', 'name', 'score', 'strand']].to_csv(
                    event_file, sep='\t', header=False, index=False, mode='a'
                )
                print(f"    {event}: {len(event_df)} events")

print("\nBED files created successfully!")
PYEOF

export BASE_DIR

# =============================================================================
# 2. Link BigWig Files (created in step 3)
# =============================================================================
echo ""
echo "=== Linking BigWig coverage files ==="

# Sample mapping
declare -A SAMPLE_NAMES
SAMPLE_NAMES[1]="EMX1_wt"
SAMPLE_NAMES[2]="EMX1_mut"
SAMPLE_NAMES[3]="Nestin_wt"
SAMPLE_NAMES[4]="Nestin_mut"

# Create symlinks to existing BigWig files
for sample_id in "${!SAMPLE_NAMES[@]}"; do
    sample_name="${SAMPLE_NAMES[$sample_id]}"
    src_bigwig="${BIGWIG_DIR}/${sample_id}_coverage.bw"
    dst_bigwig="${OUTPUT_DIR}/bigwig/${sample_name}.bw"

    if [[ -f "${src_bigwig}" ]] && [[ ! -f "${dst_bigwig}" ]]; then
        ln -sf "${src_bigwig}" "${dst_bigwig}"
        echo "  Linked: ${sample_name}.bw"
    elif [[ -f "${dst_bigwig}" ]]; then
        echo "  Already exists: ${sample_name}.bw"
    else
        echo "  Warning: Source not found: ${src_bigwig}"
    fi
done

# =============================================================================
# 3. Create IGV Session XML File
# =============================================================================
echo ""
echo "=== Creating IGV session file ==="

cat > ${OUTPUT_DIR}/splicing_analysis.xml << 'XMLEOF'
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="mm39" hasGeneTrack="true" hasSequenceTrack="true" locus="All" version="8">
    <Resources>
XMLEOF

# Add BigWig resources
for sample_name in EMX1_wt EMX1_mut Nestin_wt Nestin_mut; do
    echo "        <Resource path=\"bigwig/${sample_name}.bw\"/>" >> ${OUTPUT_DIR}/splicing_analysis.xml
done

# Add BAM resources (for sashimi plots)
for sample_id in 1 2 3 4; do
    echo "        <Resource path=\"${BAM_DIR}/${sample_id}/${sample_id}_Aligned.sortedByCoord.out.bam\"/>" >> ${OUTPUT_DIR}/splicing_analysis.xml
done

# Add BED resources (one per comparison)
for comp in EMX1_wt_vs_mut Nestin_wt_vs_mut WT_EMX1_vs_Nestin MUT_EMX1_vs_Nestin; do
    if [[ -f "${OUTPUT_DIR}/bed/${comp}_significant_splicing.bed" ]]; then
        echo "        <Resource path=\"bed/${comp}_significant_splicing.bed\"/>" >> ${OUTPUT_DIR}/splicing_analysis.xml
    fi
done

cat >> ${OUTPUT_DIR}/splicing_analysis.xml << 'XMLEOF'
    </Resources>
    <Panel name="DataPanel" height="400">
        <!-- Coverage tracks will be loaded here -->
    </Panel>
    <Panel name="FeaturePanel" height="150">
        <!-- BED tracks will be loaded here -->
    </Panel>
    <PanelLayout dividerFractions="0.72"/>
</Session>
XMLEOF

echo "  Created: splicing_analysis.xml"

# =============================================================================
# 4. Create Top Splicing Regions Navigation File
# =============================================================================
echo ""
echo "=== Creating top splicing regions file ==="

python3 << 'PYEOF2'
import os
import pandas as pd
from pathlib import Path

base_dir = os.environ.get('BASE_DIR', '/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649')
splicing_dir = Path(base_dir) / 'results/05_splicing'
output_dir = Path(base_dir) / 'results/08_igv'

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

# Collect top events across all comparisons
all_top_events = []

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

            # Filter and sort by |dPSI|
            sig = df[
                (abs(df['IncLevelDifference']) > DPSI_THRESHOLD) &
                (df['min_total_reads'] >= MIN_READS)
            ].copy()
            sig['abs_dpsi'] = abs(sig['IncLevelDifference'])

            # Get top 5 by absolute dPSI for each event type
            top = sig.nlargest(5, 'abs_dpsi')

            for _, row in top.iterrows():
                chrom = row['chr']

                # Get region with padding for visualization
                if event == 'SE':
                    start = max(1, int(row.get('upstreamES', row['exonStart_0base'])) - 500)
                    end = int(row.get('downstreamEE', row['exonEnd'])) + 500
                elif event == 'RI':
                    start = max(1, int(row['riExonStart_0base']) - 500)
                    end = int(row['riExonEnd']) + 500
                elif event in ['A5SS', 'A3SS']:
                    start = max(1, int(row.get('longExonStart_0base', row.get('shortES', 0))) - 500)
                    end = int(row.get('longExonEnd', row.get('shortEE', 0))) + 500
                elif event == 'MXE':
                    start = max(1, int(row.get('upstreamES', row['1stExonStart_0base'])) - 500)
                    end = int(row.get('downstreamEE', row['2ndExonEnd'])) + 500

                gene = str(row.get('geneSymbol', row.get('GeneID', 'Unknown'))).replace('"', '')
                dpsi = row['IncLevelDifference']

                all_top_events.append({
                    'region': f"{chrom}:{start}-{end}",
                    'description': f"{gene}_{event}_{comp}_dPSI={dpsi:.2f}",
                    'abs_dpsi': abs(dpsi),
                    'gene': gene,
                    'event': event,
                    'comparison': comp
                })
        except Exception as e:
            pass

# Sort by absolute dPSI and take top 100
all_top_events = sorted(all_top_events, key=lambda x: -x['abs_dpsi'])[:100]

# Write regions file (IGV batch format)
with open(output_dir / 'top_splicing_regions.txt', 'w') as f:
    f.write("# Top differential splicing regions for IGV\n")
    f.write("# Use: View > Go to regions in file\n")
    f.write("# Format: region<tab>description\n")
    f.write("#\n")
    for evt in all_top_events:
        f.write(f"{evt['region']}\t{evt['description']}\n")

# Write as BED for region navigation
with open(output_dir / 'bed/top_regions.bed', 'w') as f:
    f.write('track name="Top_Splicing_Regions" description="Top differential splicing events by |dPSI|"\n')
    for evt in all_top_events:
        region = evt['region']
        chrom, coords = region.split(':')
        start, end = coords.split('-')
        f.write(f"{chrom}\t{start}\t{end}\t{evt['description']}\n")

print(f"Created regions file with {len(all_top_events)} top events")
PYEOF2

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "============================================"
echo "IGV Preparation Complete!"
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "Files created:"
echo "  - bed/                      : Color-coded BED files for splicing events"
echo "  - bed/*_SE.bed              : Skipped Exon events (Red)"
echo "  - bed/*_A5SS.bed            : Alt 5' Splice Site events (Green)"
echo "  - bed/*_A3SS.bed            : Alt 3' Splice Site events (Blue)"
echo "  - bed/*_MXE.bed             : Mutually Exclusive Exon events (Orange)"
echo "  - bed/*_RI.bed              : Retained Intron events (Purple)"
echo "  - bed/top_regions.bed       : Top splicing regions"
echo "  - bigwig/                   : Coverage tracks (symlinks)"
echo "  - splicing_analysis.xml     : IGV session file"
echo "  - top_splicing_regions.txt  : Quick navigation regions"
echo ""
echo "To use in IGV:"
echo "  1. Open IGV and set genome to mm39 (GRCm39)"
echo "  2. File > Load Session > splicing_analysis.xml"
echo "  3. Or manually load BED/BigWig files"
echo "  4. Right-click on BAM track > Sashimi Plot"
echo ""
echo "Color coding:"
echo "  SE   = Red (255,0,0)      - Skipped Exon"
echo "  A5SS = Green (0,255,0)    - Alternative 5' Splice Site"
echo "  A3SS = Blue (0,0,255)     - Alternative 3' Splice Site"
echo "  MXE  = Orange (255,165,0) - Mutually Exclusive Exons"
echo "  RI   = Purple (128,0,128) - Retained Intron"
echo ""
echo "End time: $(date)"
echo "============================================"
