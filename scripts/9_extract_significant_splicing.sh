#!/bin/bash
# Extract significant differential splicing events (FDR < 0.05) from rMATS output
# Creates IGV-friendly CSV files for easy examination

set -euo pipefail

# Configuration
INPUT_DIR="results/05_splicing"
OUTPUT_DIR="results/05_splicing/significant_events"
FDR_THRESHOLD=0.05

# Comparisons and event types
COMPARISONS=("EMX1_wt_vs_mut" "Nestin_wt_vs_mut" "WT_EMX1_vs_Nestin" "MUT_EMX1_vs_Nestin")
EVENT_TYPES=("SE" "A5SS" "A3SS" "RI" "MXE")

# Create output directory
mkdir -p "$OUTPUT_DIR"

# CSV header
HEADER="comparison,event_type,gene_id,gene_symbol,chr,strand,coord1_start,coord1_end,coord2_start,coord2_end,coord3_start,coord3_end,pvalue,FDR,IncLevel_Sample1,IncLevel_Sample2,IncLevelDifference,IGV_location"

# Initialize master file
echo "$HEADER" > "$OUTPUT_DIR/all_significant_events.csv"

# Initialize summary file
echo "Significant Splicing Events Summary (FDR < $FDR_THRESHOLD)" > "$OUTPUT_DIR/significant_events_summary.txt"
echo "Generated: $(date)" >> "$OUTPUT_DIR/significant_events_summary.txt"
echo "============================================================" >> "$OUTPUT_DIR/significant_events_summary.txt"
echo "" >> "$OUTPUT_DIR/significant_events_summary.txt"

# Function to extract significant events from a single file
extract_significant() {
    local comparison=$1
    local event_type=$2
    local input_file=$3

    if [[ ! -f "$input_file" ]]; then
        echo "Warning: $input_file not found, skipping..."
        return
    fi

    # Process file: filter by FDR < threshold and format as CSV
    # Columns: 2=GeneID, 3=geneSymbol, 4=chr, 5=strand, 6-11=coords, 19=PValue, 20=FDR, 21=IncLevel1, 22=IncLevel2, 23=IncLevelDiff
    awk -F'\t' -v comp="$comparison" -v etype="$event_type" -v thresh="$FDR_THRESHOLD" '
    NR > 1 && $20 != "NA" && $20 + 0 < thresh {
        # Create IGV location (using first coordinate pair for simplicity)
        # Use 1-based coordinates for IGV
        igv_loc = $4 ":" ($6+1) "-" $7

        # Print CSV line
        printf "%s,%s,%s,%s,%s,%s,%d,%d,%d,%d,%d,%d,%s,%s,%s,%s,%s,%s\n",
            comp, etype, $2, $3, $4, $5,
            $6+1, $7, $8+1, $9, $10+1, $11,
            $19, $20, $21, $22, $23, igv_loc
    }' "$input_file"
}

# Process each comparison
for comparison in "${COMPARISONS[@]}"; do
    echo "Processing $comparison..."

    # Initialize per-comparison file
    comp_file="$OUTPUT_DIR/${comparison}_significant.csv"
    echo "$HEADER" > "$comp_file"

    comp_total=0

    for event_type in "${EVENT_TYPES[@]}"; do
        input_file="$INPUT_DIR/$comparison/${event_type}.MATS.JCEC.txt"

        if [[ -f "$input_file" ]]; then
            # Extract significant events
            events=$(extract_significant "$comparison" "$event_type" "$input_file")

            if [[ -n "$events" ]]; then
                # Count events
                event_count=$(echo "$events" | wc -l)
                comp_total=$((comp_total + event_count))

                # Add to comparison file (sorted by abs(IncLevelDifference))
                echo "$events" >> "$comp_file"

                # Add to master file
                echo "$events" >> "$OUTPUT_DIR/all_significant_events.csv"

                echo "  $event_type: $event_count significant events"
            else
                echo "  $event_type: 0 significant events"
            fi
        else
            echo "  $event_type: file not found"
        fi
    done

    # Sort comparison file by absolute IncLevelDifference (column 17)
    if [[ -f "$comp_file" ]]; then
        # Keep header, sort rest by abs(IncLevelDiff) descending
        (head -1 "$comp_file" && tail -n +2 "$comp_file" | awk -F',' '{
            val = $17
            if (val < 0) val = -val
            print val "\t" $0
        }' | sort -t$'\t' -k1 -rn | cut -f2-) > "${comp_file}.tmp"
        mv "${comp_file}.tmp" "$comp_file"
    fi

    echo "  Total: $comp_total significant events for $comparison"
    echo "" >> "$OUTPUT_DIR/significant_events_summary.txt"
    echo "$comparison: $comp_total total significant events" >> "$OUTPUT_DIR/significant_events_summary.txt"
done

# Sort master file by comparison, then by abs(IncLevelDifference)
echo ""
echo "Sorting master file..."
(head -1 "$OUTPUT_DIR/all_significant_events.csv" && tail -n +2 "$OUTPUT_DIR/all_significant_events.csv" | sort -t',' -k1,1 -k17,17rn) > "$OUTPUT_DIR/all_significant_events.csv.tmp"
mv "$OUTPUT_DIR/all_significant_events.csv.tmp" "$OUTPUT_DIR/all_significant_events.csv"

# Final summary
total_events=$(tail -n +2 "$OUTPUT_DIR/all_significant_events.csv" | wc -l)
echo ""
echo "============================================================" >> "$OUTPUT_DIR/significant_events_summary.txt"
echo "" >> "$OUTPUT_DIR/significant_events_summary.txt"
echo "Total significant events across all comparisons: $total_events" >> "$OUTPUT_DIR/significant_events_summary.txt"
echo "" >> "$OUTPUT_DIR/significant_events_summary.txt"

# Add event type breakdown
echo "Event type breakdown:" >> "$OUTPUT_DIR/significant_events_summary.txt"
for event_type in "${EVENT_TYPES[@]}"; do
    count=$(tail -n +2 "$OUTPUT_DIR/all_significant_events.csv" | awk -F',' -v et="$event_type" '$2 == et' | wc -l)
    echo "  $event_type: $count" >> "$OUTPUT_DIR/significant_events_summary.txt"
done

echo ""
echo "Done! Output files:"
echo "  - Per-comparison files: $OUTPUT_DIR/{comparison}_significant.csv"
echo "  - Master file: $OUTPUT_DIR/all_significant_events.csv"
echo "  - Summary: $OUTPUT_DIR/significant_events_summary.txt"
echo ""
echo "Total significant events: $total_events"
