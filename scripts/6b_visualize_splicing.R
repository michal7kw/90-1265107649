#!/usr/bin/env Rscript

#===============================================================================
# SCRIPT: 6b_visualize_splicing.R
# PURPOSE: Visualize rMATS Differential Splicing Results
#
# DESCRIPTION:
# Creates summary plots and tables from rMATS output files.
#
# IMPORTANT - MONOPLICATE DESIGN LIMITATIONS:
# This experiment has n=1 per condition (no biological replicates).
# Splicing results should be treated as EXPLORATORY only.
# Filtering criteria for this analysis:
#   - |deltaPSI| > 0.1 (10% difference in inclusion)
#   - Minimum total junction reads > 20 (for reliability)
# Even with these filters, false positives are expected.
# Validate candidate events experimentally (RT-PCR, etc.)
#
# USAGE:
# Rscript 6b_visualize_splicing.R
#===============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
})

cat("=== Visualizing Splicing Results ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
SPLICING_DIR <- file.path(BASE_DIR, "results/05_splicing")
OUTPUT_DIR <- file.path(SPLICING_DIR, "visualizations")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Comparisons to analyze
comparisons <- c("EMX1_wt_vs_mut", "Nestin_wt_vs_mut", "WT_EMX1_vs_Nestin", "MUT_EMX1_vs_Nestin")
event_types <- c("SE", "A5SS", "A3SS", "MXE", "RI")

# Filtering thresholds (stricter for monoplicate data)
PSI_THRESHOLD <- 0.1  # Minimum |deltaPSI|
MIN_READS <- 20       # Minimum total junction reads for reliability

#-------------------------------------------------------------------------------
# Function to read and process rMATS output
#-------------------------------------------------------------------------------
read_rmats_results <- function(comparison, event_type) {
    file_path <- file.path(SPLICING_DIR, comparison, paste0(event_type, ".MATS.JC.txt"))

    if (!file.exists(file_path)) {
        return(NULL)
    }

    df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # Add metadata
    df$comparison <- comparison
    df$event_type <- event_type

    # Calculate total junction reads for filtering
    # IJC = Inclusion Junction Count, SJC = Skipping Junction Count
    df$total_reads_s1 <- df$IJC_SAMPLE_1 + df$SJC_SAMPLE_1
    df$total_reads_s2 <- df$IJC_SAMPLE_2 + df$SJC_SAMPLE_2
    df$min_total_reads <- pmin(df$total_reads_s1, df$total_reads_s2)

    # Filter significant events (stricter for monoplicate data)
    # Require both PSI difference AND minimum read support
    df$significant <- abs(df$IncLevelDifference) > PSI_THRESHOLD &
                      df$min_total_reads >= MIN_READS

    return(df)
}

#-------------------------------------------------------------------------------
# Collect all results
#-------------------------------------------------------------------------------
cat("Loading rMATS results...\n")

all_results <- list()

for (comp in comparisons) {
    for (event in event_types) {
        result <- read_rmats_results(comp, event)
        if (!is.null(result)) {
            all_results[[paste(comp, event, sep = "_")]] <- result
        }
    }
}

if (length(all_results) == 0) {
    cat("No rMATS results found. Please run differential splicing analysis first.\n")
    quit(status = 0)
}

# Combine all results
combined_df <- bind_rows(all_results)
cat("Total events loaded:", nrow(combined_df), "\n")
cat("Significant events (|deltaPSI| >", PSI_THRESHOLD, "):", sum(combined_df$significant), "\n\n")

#-------------------------------------------------------------------------------
# 1. Summary bar plot
#-------------------------------------------------------------------------------
cat("Creating summary plots...\n")

summary_df <- combined_df %>%
    filter(significant) %>%
    group_by(comparison, event_type) %>%
    summarise(count = n(), .groups = "drop")

# Add zeros for missing combinations
full_grid <- expand.grid(comparison = comparisons, event_type = event_types)
summary_df <- full_grid %>%
    left_join(summary_df, by = c("comparison", "event_type")) %>%
    replace_na(list(count = 0))

p1 <- ggplot(summary_df, aes(x = comparison, y = count, fill = event_type)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0("Significant Splicing Events"),
         subtitle = paste0("Filters: |deltaPSI| > ", PSI_THRESHOLD,
                          " AND min junction reads >= ", MIN_READS,
                          "\n(Monoplicate data - interpret with caution)"),
         x = "Comparison",
         y = "Number of Events",
         fill = "Event Type") +
    scale_fill_brewer(palette = "Set2")

ggsave(file.path(OUTPUT_DIR, "splicing_events_summary.pdf"), p1, width = 10, height = 6)
ggsave(file.path(OUTPUT_DIR, "splicing_events_summary.png"), p1, width = 10, height = 6, dpi = 300)

#-------------------------------------------------------------------------------
# 2. Volcano-like plot (deltaPSI distribution)
#-------------------------------------------------------------------------------
p2 <- ggplot(combined_df, aes(x = IncLevelDifference, fill = event_type)) +
    geom_histogram(bins = 50, alpha = 0.7) +
    facet_wrap(~ comparison, ncol = 2) +
    geom_vline(xintercept = c(-PSI_THRESHOLD, PSI_THRESHOLD), linetype = "dashed", color = "red") +
    theme_bw() +
    labs(title = "Distribution of PSI Differences",
         x = "Delta PSI (Inclusion Level Difference)",
         y = "Count",
         fill = "Event Type") +
    scale_fill_brewer(palette = "Set2")

ggsave(file.path(OUTPUT_DIR, "deltaPSI_distribution.pdf"), p2, width = 12, height = 10)
ggsave(file.path(OUTPUT_DIR, "deltaPSI_distribution.png"), p2, width = 12, height = 10, dpi = 300)

#-------------------------------------------------------------------------------
# 3. Event type breakdown
#-------------------------------------------------------------------------------
event_summary <- combined_df %>%
    group_by(event_type) %>%
    summarise(
        total = n(),
        significant = sum(significant),
        pct_significant = round(100 * significant / total, 1)
    )

cat("\nEvent Type Summary:\n")
print(event_summary)

write.csv(event_summary, file.path(OUTPUT_DIR, "event_type_summary.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# 4. Top significant events table
#-------------------------------------------------------------------------------
sig_events <- combined_df %>%
    filter(significant) %>%
    arrange(desc(abs(IncLevelDifference))) %>%
    select(comparison, event_type, GeneID, geneSymbol, chr, strand,
           IncLevelDifference, min_total_reads,
           IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_2) %>%
    head(100)

write.csv(sig_events, file.path(OUTPUT_DIR, "top_significant_events.csv"), row.names = FALSE)

cat("\nFiltering criteria applied:\n")
cat("  - |deltaPSI| >", PSI_THRESHOLD, "\n")
cat("  - Minimum junction reads >=", MIN_READS, "\n")
cat("Top significant events saved:", nrow(sig_events), "\n")

#-------------------------------------------------------------------------------
# 5. Gene-level summary
#-------------------------------------------------------------------------------
gene_summary <- combined_df %>%
    filter(significant) %>%
    group_by(geneSymbol, comparison) %>%
    summarise(
        n_events = n(),
        event_types = paste(unique(event_type), collapse = ","),
        max_deltaPSI = max(abs(IncLevelDifference)),
        .groups = "drop"
    ) %>%
    arrange(desc(n_events))

write.csv(gene_summary, file.path(OUTPUT_DIR, "gene_level_summary.csv"), row.names = FALSE)

cat("\nTop genes with multiple splicing events:\n")
print(head(gene_summary, 20))

#-------------------------------------------------------------------------------
# 6. Overlap analysis between comparisons
#-------------------------------------------------------------------------------
cat("\nAnalyzing overlaps between comparisons...\n")

# Get significant genes for each comparison
sig_genes_by_comp <- combined_df %>%
    filter(significant) %>%
    group_by(comparison) %>%
    summarise(genes = list(unique(geneSymbol)), .groups = "drop")

# Create overlap matrix
overlap_matrix <- matrix(0, nrow = length(comparisons), ncol = length(comparisons),
                         dimnames = list(comparisons, comparisons))

for (i in 1:length(comparisons)) {
    for (j in 1:length(comparisons)) {
        genes_i <- combined_df %>%
            filter(comparison == comparisons[i], significant) %>%
            pull(geneSymbol) %>% unique()
        genes_j <- combined_df %>%
            filter(comparison == comparisons[j], significant) %>%
            pull(geneSymbol) %>% unique()
        overlap_matrix[i, j] <- length(intersect(genes_i, genes_j))
    }
}

pdf(file.path(OUTPUT_DIR, "comparison_overlap_heatmap.pdf"), width = 8, height = 7)
pheatmap(overlap_matrix,
         display_numbers = TRUE,
         number_format = "%d",
         main = "Overlap of Genes with Significant Splicing Events",
         color = colorRampPalette(brewer.pal(9, "Blues"))(100))
dev.off()

#-------------------------------------------------------------------------------
# 7. Summary report
#-------------------------------------------------------------------------------
summary_report <- data.frame(
    Comparison = comparisons,
    Total_Events = sapply(comparisons, function(x) sum(combined_df$comparison == x)),
    Significant_Events = sapply(comparisons, function(x) sum(combined_df$comparison == x & combined_df$significant)),
    Affected_Genes = sapply(comparisons, function(x) {
        combined_df %>%
            filter(comparison == x, significant) %>%
            pull(geneSymbol) %>% unique() %>% length()
    })
)

write.csv(summary_report, file.path(OUTPUT_DIR, "splicing_analysis_summary.csv"), row.names = FALSE)

cat("\n=== Summary Report ===\n")
print(summary_report)

cat("\n=== IMPORTANT CAVEATS ===\n")
cat("This analysis was performed on MONOPLICATE data (n=1 per condition).\n")
cat("Results should be interpreted as EXPLORATORY only.\n")
cat("False positive rate is expected to be HIGH.\n")
cat("Key findings should be validated experimentally (RT-PCR, etc.)\n")

cat("\n=== Visualization Complete ===\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
