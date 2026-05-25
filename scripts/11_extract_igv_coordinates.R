#!/usr/bin/env Rscript

#===============================================================================
# SCRIPT: 11_extract_igv_coordinates.R
# PURPOSE: Extract and consolidate significant splicing events with IGV coordinates
#
# DESCRIPTION:
# Creates multiple output files for easy IGV navigation:
#   1. Master CSV with all significant events + IGV coordinates
#   2. Per-comparison navigation files
#   3. BED files per event type for IGV
#   4. Quick-copy tab-separated files for easy paste into IGV
#   5. Top hits file (high |dPSI| events)
#
# IMPORTANT - MONOPLICATE DESIGN:
# This experiment has n=1 per condition. We use:
#   - |deltaPSI| > 0.1 (not FDR)
#   - Minimum junction reads >= 20
#
# USAGE:
# Rscript 11_extract_igv_coordinates.R
#===============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})

#===============================================================================
# Configuration
#===============================================================================

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
SPLICING_DIR <- file.path(BASE_DIR, "results/05_splicing")
OUTPUT_DIR <- file.path(BASE_DIR, "results/08_igv/coordinates")

# Significance thresholds (for monoplicate data)
DPSI_THRESHOLD <- 0.1
MIN_READS <- 20

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Define comparisons and event types
COMPARISONS <- c(
    "EMX1_wt_vs_mut",
    "Nestin_wt_vs_mut",
    "WT_EMX1_vs_Nestin",
    "MUT_EMX1_vs_Nestin"
)

EVENT_TYPES <- list(
    SE = list(
        file = "SE.MATS.JC.txt",
        desc = "Skipped Exon",
        color = "255,0,0"
    ),
    A3SS = list(
        file = "A3SS.MATS.JC.txt",
        desc = "Alternative 3' Splice Site",
        color = "0,0,255"
    ),
    A5SS = list(
        file = "A5SS.MATS.JC.txt",
        desc = "Alternative 5' Splice Site",
        color = "0,255,0"
    ),
    MXE = list(
        file = "MXE.MATS.JC.txt",
        desc = "Mutually Exclusive Exons",
        color = "255,165,0"
    ),
    RI = list(
        file = "RI.MATS.JC.txt",
        desc = "Retained Intron",
        color = "128,0,128"
    )
)

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Extracting Significant Splicing Events with IGV Coordinates\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("|dPSI| threshold: %s\n", DPSI_THRESHOLD))
cat(sprintf("Minimum reads threshold: %s\n", MIN_READS))
cat("\n")

#===============================================================================
# Helper Functions
#===============================================================================

# Sum counts from comma-separated string
sum_counts <- function(x) {
    if (is.null(x) || is.na(x) || x == "") return(0)
    sum(as.numeric(unlist(strsplit(as.character(x), ","))), na.rm = TRUE)
}

# Read and filter rMATS output
read_rmats_significant <- function(comparison, event_type, event_info) {
    file_path <- file.path(SPLICING_DIR, comparison, event_info$file)

    if (!file.exists(file_path)) {
        message(sprintf("  File not found: %s", file_path))
        return(NULL)
    }

    # Read rMATS output
    df <- read.delim(file_path, stringsAsFactors = FALSE)

    # Calculate junction read counts
    if (all(c("IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2") %in% colnames(df))) {
        df$total_reads_1 <- sapply(df$IJC_SAMPLE_1, sum_counts) + sapply(df$SJC_SAMPLE_1, sum_counts)
        df$total_reads_2 <- sapply(df$IJC_SAMPLE_2, sum_counts) + sapply(df$SJC_SAMPLE_2, sum_counts)
        df$min_total_reads <- pmin(df$total_reads_1, df$total_reads_2)
    } else {
        df$min_total_reads <- MIN_READS
    }

    # Filter for significant events (using deltaPSI, not FDR)
    df_sig <- df %>%
        filter(
            abs(IncLevelDifference) >= DPSI_THRESHOLD,
            min_total_reads >= MIN_READS
        )

    if (nrow(df_sig) == 0) {
        return(NULL)
    }

    # Add metadata columns
    df_sig <- df_sig %>%
        mutate(
            Comparison = comparison,
            EventType = event_type,
            EventDescription = event_info$desc,
            DeltaPSI = IncLevelDifference,
            PSI_Group1 = IncLevel1,
            PSI_Group2 = IncLevel2,
            Direction = ifelse(IncLevelDifference > 0, "Inclusion", "Exclusion")
        )

    # Standardize column names
    df_sig <- df_sig %>%
        rename(
            GeneSymbol = geneSymbol,
            Chromosome = chr
        )

    return(df_sig)
}

# Add IGV coordinates based on event type
add_igv_coords <- function(df, event_type) {
    if (event_type == "SE") {
        df <- df %>%
            mutate(
                event_start = exonStart_0base,
                event_end = exonEnd,
                ExonLength = exonEnd - exonStart_0base
            )
    } else if (event_type %in% c("A3SS", "A5SS")) {
        df <- df %>%
            mutate(
                event_start = longExonStart_0base,
                event_end = longExonEnd,
                ExonLength = NA_real_
            )
    } else if (event_type == "MXE") {
        # Handle MXE columns (may have X prefix due to numeric names)
        if ("X1stExonStart_0base" %in% colnames(df)) {
            df <- df %>%
                mutate(
                    event_start = X1stExonStart_0base,
                    event_end = X2ndExonEnd,
                    ExonLength = NA_real_
                )
        } else if ("1stExonStart_0base" %in% colnames(df)) {
            df <- df %>%
                mutate(
                    event_start = `1stExonStart_0base`,
                    event_end = `2ndExonEnd`,
                    ExonLength = NA_real_
                )
        }
    } else if (event_type == "RI") {
        df <- df %>%
            mutate(
                event_start = riExonStart_0base,
                event_end = riExonEnd,
                ExonLength = NA_real_
            )
    }

    # Create IGV coordinate strings
    df %>%
        mutate(
            IGV_coordinate = paste0(
                Chromosome, ":",
                pmax(0, event_start - 500), "-",
                event_end + 500
            ),
            IGV_narrow = paste0(
                Chromosome, ":",
                pmax(0, event_start - 100), "-",
                event_end + 100
            )
        )
}

#===============================================================================
# Main Processing
#===============================================================================

# Collect all significant events
all_events <- list()
comparison_events <- list()

for (comparison in COMPARISONS) {
    cat(sprintf("Processing %s...\n", comparison))
    comp_events <- list()

    for (event_type in names(EVENT_TYPES)) {
        event_info <- EVENT_TYPES[[event_type]]
        df <- read_rmats_significant(comparison, event_type, event_info)

        if (!is.null(df) && nrow(df) > 0) {
            cat(sprintf("  %s: %d significant events\n", event_type, nrow(df)))
            comp_events[[event_type]] <- df
            all_events[[paste(comparison, event_type, sep = "_")]] <- df
        }
    }

    # Combine events for this comparison
    if (length(comp_events) > 0) {
        comparison_events[[comparison]] <- bind_rows(comp_events)
    }
}

#===============================================================================
# Create Output Files
#===============================================================================

cat("\nCreating output files...\n")

# Combine all events
all_events_df <- bind_rows(all_events) %>%
    arrange(Comparison, EventType, desc(abs(DeltaPSI)))

# Process each event type separately for IGV coordinates
get_all_events_with_coords <- function(df) {
    event_types <- unique(df$EventType)
    result_list <- list()

    for (et in event_types) {
        et_df <- df %>% filter(EventType == et)
        if (nrow(et_df) > 0) {
            result_list[[et]] <- add_igv_coords(et_df, et)
        }
    }

    bind_rows(result_list)
}

all_events_igv <- get_all_events_with_coords(all_events_df)

#-------------------------------------------------------------------------------
# 1. Master CSV with all significant events
#-------------------------------------------------------------------------------

master_df <- all_events_igv %>%
    select(
        Comparison, EventType, EventDescription, GeneSymbol, GeneID,
        Chromosome, strand, DeltaPSI, Direction, min_total_reads,
        event_start, event_end, ExonLength, IGV_coordinate, IGV_narrow,
        PSI_Group1, PSI_Group2
    ) %>%
    arrange(Comparison, desc(abs(DeltaPSI)))

write.csv(master_df,
          file.path(OUTPUT_DIR, "all_significant_events_IGV.csv"),
          row.names = FALSE)
cat(sprintf("  all_significant_events_IGV.csv: %d events\n", nrow(master_df)))

#-------------------------------------------------------------------------------
# 2. Simplified version (key columns only)
#-------------------------------------------------------------------------------

simple_df <- master_df %>%
    select(
        Gene = GeneSymbol,
        Comparison, EventType, DeltaPSI, Direction,
        Chromosome, Start = event_start, End = event_end,
        IGV_coordinate
    ) %>%
    mutate(DeltaPSI = round(DeltaPSI, 3))

write.csv(simple_df,
          file.path(OUTPUT_DIR, "significant_events_simple.csv"),
          row.names = FALSE)
cat(sprintf("  significant_events_simple.csv: %d events\n", nrow(simple_df)))

#-------------------------------------------------------------------------------
# 3. Per-comparison files
#-------------------------------------------------------------------------------

for (comparison in names(comparison_events)) {
    df <- comparison_events[[comparison]]
    df_igv <- get_all_events_with_coords(df)

    filename <- sprintf("%s_significant_events.csv", comparison)
    write.csv(df_igv, file.path(OUTPUT_DIR, filename), row.names = FALSE)
    cat(sprintf("  %s: %d events\n", filename, nrow(df_igv)))

    # Also create simple navigation file (tab-separated for copy-paste)
    nav_df <- df_igv %>%
        mutate(
            Label = paste0(GeneSymbol, " | ", EventType, " | dPSI=", round(DeltaPSI, 2))
        ) %>%
        select(Label, IGV_coordinate) %>%
        distinct()

    nav_file <- sprintf("IGV_nav_%s.txt", comparison)
    write.table(nav_df, file.path(OUTPUT_DIR, nav_file),
                sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}
cat("  Per-comparison IGV navigation files created\n")

#-------------------------------------------------------------------------------
# 4. Gene-level summary
#-------------------------------------------------------------------------------

genes_summary <- master_df %>%
    group_by(GeneSymbol) %>%
    summarise(
        Comparisons = paste(unique(Comparison), collapse = ","),
        EventTypes = paste(unique(EventType), collapse = ","),
        TotalEvents = n(),
        Mean_dPSI = round(mean(DeltaPSI), 3),
        Max_abs_dPSI = round(max(abs(DeltaPSI)), 3),
        .groups = "drop"
    ) %>%
    arrange(desc(TotalEvents))

write.csv(genes_summary,
          file.path(OUTPUT_DIR, "genes_with_splicing_events.csv"),
          row.names = FALSE)
cat(sprintf("  genes_with_splicing_events.csv: %d genes\n", nrow(genes_summary)))

#-------------------------------------------------------------------------------
# 5. Summary by comparison and event type
#-------------------------------------------------------------------------------

summary_stats <- master_df %>%
    group_by(Comparison, EventType) %>%
    summarise(
        Count = n(),
        Mean_dPSI = round(mean(DeltaPSI), 3),
        Mean_abs_dPSI = round(mean(abs(DeltaPSI)), 3),
        Inclusion_count = sum(Direction == "Inclusion"),
        Exclusion_count = sum(Direction == "Exclusion"),
        .groups = "drop"
    )

write.csv(summary_stats,
          file.path(OUTPUT_DIR, "summary_by_comparison_eventtype.csv"),
          row.names = FALSE)
cat(sprintf("  summary_by_comparison_eventtype.csv: %d rows\n", nrow(summary_stats)))

#-------------------------------------------------------------------------------
# 6. Quick-copy navigation file (all events)
#-------------------------------------------------------------------------------

quick_nav <- master_df %>%
    mutate(
        Label = paste0(GeneSymbol, " | ", EventType, " | dPSI=", round(DeltaPSI, 2),
                      " | ", Comparison)
    ) %>%
    select(Label, IGV_coordinate) %>%
    distinct()

write.table(quick_nav, file.path(OUTPUT_DIR, "IGV_coordinates_all.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
cat(sprintf("  IGV_coordinates_all.txt: %d coordinates\n", nrow(quick_nav)))

#-------------------------------------------------------------------------------
# 7. Top hits file (|dPSI| >= 0.2)
#-------------------------------------------------------------------------------

top_hits <- master_df %>%
    filter(abs(DeltaPSI) >= 0.2) %>%
    arrange(desc(abs(DeltaPSI)))

write.csv(top_hits,
          file.path(OUTPUT_DIR, "IGV_top_hits_dPSI_0.2.csv"),
          row.names = FALSE)
cat(sprintf("  IGV_top_hits_dPSI_0.2.csv: %d events with |dPSI| >= 0.2\n", nrow(top_hits)))

# Quick nav for top hits
if (nrow(top_hits) > 0) {
    top_nav <- top_hits %>%
        mutate(
            Label = paste0(GeneSymbol, " | ", EventType, " | dPSI=", round(DeltaPSI, 2))
        ) %>%
        select(Label, IGV_coordinate) %>%
        distinct()

    write.table(top_nav, file.path(OUTPUT_DIR, "IGV_top_hits.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}

#-------------------------------------------------------------------------------
# 8. BED files per event type (for IGV color-coded loading)
#-------------------------------------------------------------------------------

bed_dir <- file.path(OUTPUT_DIR, "bed")
dir.create(bed_dir, showWarnings = FALSE)

for (event_type in names(EVENT_TYPES)) {
    event_df <- master_df %>%
        filter(EventType == event_type)

    if (nrow(event_df) > 0) {
        bed_file <- file.path(bed_dir, paste0("all_", event_type, ".bed"))
        color <- EVENT_TYPES[[event_type]]$color

        # BED6 format
        bed_data <- event_df %>%
            mutate(
                score = pmin(1000, round(abs(DeltaPSI) * 1000)),
                name = paste0(GeneSymbol, "_", Comparison, "_dPSI=", round(DeltaPSI, 2))
            ) %>%
            select(
                chrom = Chromosome,
                start = event_start,
                end = event_end,
                name,
                score,
                strand
            )

        # Write with track header
        con <- file(bed_file, "w")
        writeLines(sprintf('track name="%s_events" description="%s events" color="%s"',
                          event_type, EVENT_TYPES[[event_type]]$desc, gsub(",", " ", color)),
                  con)
        close(con)

        write.table(bed_data, bed_file, sep = "\t", row.names = FALSE,
                    col.names = FALSE, quote = FALSE, append = TRUE)

        cat(sprintf("  bed/all_%s.bed: %d events\n", event_type, nrow(event_df)))
    }
}

#===============================================================================
# Print Summary
#===============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat("\nTotal significant events by comparison:\n")
master_df %>%
    group_by(Comparison) %>%
    summarise(Count = n(), .groups = "drop") %>%
    arrange(desc(Count)) %>%
    print()

cat("\nTotal significant events by event type:\n")
master_df %>%
    group_by(EventType) %>%
    summarise(Count = n(), .groups = "drop") %>%
    arrange(desc(Count)) %>%
    print()

cat("\nTop 10 genes with most splicing events:\n")
genes_summary %>%
    head(10) %>%
    select(GeneSymbol, TotalEvents, EventTypes) %>%
    print()

cat("\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("\nFiles created:\n")
cat("  - all_significant_events_IGV.csv : Complete data with IGV coordinates\n")
cat("  - significant_events_simple.csv  : Simplified key columns\n")
cat("  - {Comparison}_significant_events.csv : Per-comparison files\n")
cat("  - IGV_nav_{Comparison}.txt       : Copy-paste navigation files\n")
cat("  - IGV_coordinates_all.txt        : All events for quick navigation\n")
cat("  - IGV_top_hits*.csv/txt          : High-confidence events (|dPSI| >= 0.2)\n")
cat("  - bed/all_{EventType}.bed        : Event-type-specific BED files\n")
cat("\nTo use in IGV:\n")
cat("  1. Open IGV and set genome to mm39\n")
cat("  2. Copy coordinate from *_coordinates*.txt files\n")
cat("  3. Paste into IGV search box\n")
cat("  4. Or load BED files for color-coded visualization\n")
cat("\nDone!\n")
