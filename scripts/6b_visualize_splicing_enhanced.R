#!/usr/bin/env Rscript

#===============================================================================
# SCRIPT: 6b_visualize_splicing_enhanced.R
# PURPOSE: Enhanced Visualization of rMATS Differential Splicing Results
#
# DESCRIPTION:
# Creates comprehensive visualizations for splicing analysis results including:
#   1. Event type distributions (stacked bars, pie charts)
#   2. deltaPSI density and distribution analysis
#   3. Directionality analysis (inclusion vs exclusion)
#   4. Venn diagram of overlapping genes across comparisons
#   5. Enhanced summary tables with IGV coordinates
#
# IMPORTANT - MONOPLICATE DESIGN LIMITATIONS:
# This experiment has n=1 per condition (no biological replicates).
# - PCA analysis is NOT performed (meaningless with n=1)
# - FDR values are NOT used for filtering (unreliable without replicates)
# - Instead, we use: |deltaPSI| > 0.1 AND minimum junction reads >= 20
# - Results should be treated as EXPLORATORY only
# - Validate candidate events experimentally (RT-PCR, etc.)
#
# USAGE:
# Rscript 6b_visualize_splicing_enhanced.R
#===============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(VennDiagram)
    library(grid)
    library(gridExtra)
    library(patchwork)
})

cat("=== Enhanced Splicing Visualization ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

#===============================================================================
# Configuration
#===============================================================================

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
SPLICING_DIR <- file.path(BASE_DIR, "results/05_splicing")
OUTPUT_DIR <- file.path(SPLICING_DIR, "visualizations_enhanced")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Comparisons for this project
comparisons <- c("EMX1_wt_vs_mut", "Nestin_wt_vs_mut", "WT_EMX1_vs_Nestin", "MUT_EMX1_vs_Nestin")
event_types <- c("SE", "A5SS", "A3SS", "MXE", "RI")

# Filtering thresholds (stricter for monoplicate data)
PSI_THRESHOLD <- 0.1   # Minimum |deltaPSI|
MIN_READS <- 20        # Minimum total junction reads for reliability

# Color palettes
event_colors <- c(
    SE = "#E41A1C",     # Red
    A5SS = "#4DAF4A",   # Green
    A3SS = "#377EB8",   # Blue
    MXE = "#FF7F00",    # Orange
    RI = "#984EA3"      # Purple
)

comparison_colors <- c(
    EMX1_wt_vs_mut = "#E41A1C",
    Nestin_wt_vs_mut = "#377EB8",
    WT_EMX1_vs_Nestin = "#4DAF4A",
    MUT_EMX1_vs_Nestin = "#984EA3"
)

#===============================================================================
# Helper Functions
#===============================================================================

# Function to sum counts from comma-separated string
sum_counts <- function(x) {
    if (is.null(x) || is.na(x) || x == "") return(0)
    sum(as.numeric(unlist(strsplit(as.character(x), ","))), na.rm = TRUE)
}

# Function to read and process rMATS output
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
    if (all(c("IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2") %in% colnames(df))) {
        df$total_reads_s1 <- sapply(df$IJC_SAMPLE_1, sum_counts) + sapply(df$SJC_SAMPLE_1, sum_counts)
        df$total_reads_s2 <- sapply(df$IJC_SAMPLE_2, sum_counts) + sapply(df$SJC_SAMPLE_2, sum_counts)
        df$min_total_reads <- pmin(df$total_reads_s1, df$total_reads_s2)
    } else {
        df$min_total_reads <- MIN_READS
    }

    # Define significance based on deltaPSI and read count (NOT FDR for monoplicate)
    df$significant <- abs(df$IncLevelDifference) > PSI_THRESHOLD &
                      df$min_total_reads >= MIN_READS

    # Add direction
    df$direction <- ifelse(df$IncLevelDifference > 0, "Inclusion", "Exclusion")

    return(df)
}

#===============================================================================
# Load All Results
#===============================================================================

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
cat("Significant events (|deltaPSI| >", PSI_THRESHOLD, ", reads >=", MIN_READS, "):",
    sum(combined_df$significant), "\n\n")

#===============================================================================
# 1. Event Type Distribution - Stacked Bar Charts
#===============================================================================

cat("Creating event distribution plots...\n")

# Summarize events
summary_df <- combined_df %>%
    group_by(comparison, event_type) %>%
    summarise(
        total = n(),
        significant = sum(significant),
        .groups = "drop"
    )

# Fill missing combinations with zeros
full_grid <- expand.grid(comparison = comparisons, event_type = event_types, stringsAsFactors = FALSE)
summary_df <- full_grid %>%
    left_join(summary_df, by = c("comparison", "event_type")) %>%
    replace_na(list(total = 0, significant = 0))

# 1a. Total events by type (stacked)
p1a <- ggplot(summary_df, aes(x = comparison, y = total, fill = event_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = event_colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
    labs(
        title = "Distribution of All Splicing Events by Type",
        subtitle = "All detected events (before filtering)",
        x = "Comparison",
        y = "Number of Events",
        fill = "Event Type"
    )

ggsave(file.path(OUTPUT_DIR, "event_distribution_total.pdf"), p1a, width = 10, height = 7)
ggsave(file.path(OUTPUT_DIR, "event_distribution_total.png"), p1a, width = 10, height = 7, dpi = 300)

# 1b. Significant events by type (stacked)
p1b <- ggplot(summary_df, aes(x = comparison, y = significant, fill = event_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = event_colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
    labs(
        title = "Distribution of Significant Splicing Events by Type",
        subtitle = paste0("Filters: |deltaPSI| > ", PSI_THRESHOLD, " AND min reads >= ", MIN_READS,
                         "\n(Monoplicate data - interpret with caution)"),
        x = "Comparison",
        y = "Number of Significant Events",
        fill = "Event Type"
    )

ggsave(file.path(OUTPUT_DIR, "event_distribution_significant.pdf"), p1b, width = 10, height = 7)
ggsave(file.path(OUTPUT_DIR, "event_distribution_significant.png"), p1b, width = 10, height = 7, dpi = 300)

# 1c. Side-by-side comparison (dodged)
p1c <- ggplot(summary_df, aes(x = comparison, y = significant, fill = event_type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = event_colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
        title = "Significant Splicing Events by Type (Side-by-Side)",
        x = "Comparison",
        y = "Number of Events",
        fill = "Event Type"
    )

ggsave(file.path(OUTPUT_DIR, "event_distribution_dodged.pdf"), p1c, width = 12, height = 7)

#===============================================================================
# 2. Pie Charts per Comparison
#===============================================================================

cat("Creating pie charts...\n")

for (comp in comparisons) {
    comp_data <- summary_df %>%
        filter(comparison == comp, significant > 0)

    if (nrow(comp_data) > 0 && sum(comp_data$significant) > 0) {
        comp_data <- comp_data %>%
            mutate(
                percentage = significant / sum(significant) * 100,
                label = paste0(event_type, "\n", significant, " (", round(percentage, 1), "%)")
            )

        p_pie <- ggplot(comp_data, aes(x = "", y = significant, fill = event_type)) +
            geom_bar(stat = "identity", width = 1) +
            coord_polar("y", start = 0) +
            scale_fill_manual(values = event_colors) +
            theme_void() +
            theme(legend.position = "right") +
            labs(
                title = paste0("Event Type Distribution: ", comp),
                subtitle = paste0("Total significant events: ", sum(comp_data$significant)),
                fill = "Event Type"
            ) +
            geom_text(aes(label = ifelse(significant > 0, significant, "")),
                     position = position_stack(vjust = 0.5), color = "white", fontface = "bold")

        ggsave(file.path(OUTPUT_DIR, paste0("pie_chart_", comp, ".pdf")), p_pie, width = 8, height = 6)
    }
}

#===============================================================================
# 3. deltaPSI Distribution Analysis
#===============================================================================

cat("Creating deltaPSI distribution plots...\n")

sig_events <- combined_df %>% filter(significant)

if (nrow(sig_events) > 0) {
    # 3a. Density plot by event type
    p3a <- ggplot(sig_events, aes(x = IncLevelDifference, color = event_type)) +
        geom_density(linewidth = 1, alpha = 0.7) +
        facet_wrap(~ comparison, ncol = 2) +
        scale_color_manual(values = event_colors) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
        geom_vline(xintercept = c(-PSI_THRESHOLD, PSI_THRESHOLD),
                   linetype = "dotted", color = "grey50") +
        theme_bw() +
        labs(
            title = "deltaPSI Distribution by Event Type",
            subtitle = "Significant events only",
            x = "deltaPSI (Inclusion Level Difference)",
            y = "Density",
            color = "Event Type"
        )

    ggsave(file.path(OUTPUT_DIR, "dpsi_density_by_type.pdf"), p3a, width = 12, height = 10)
    ggsave(file.path(OUTPUT_DIR, "dpsi_density_by_type.png"), p3a, width = 12, height = 10, dpi = 300)

    # 3b. Box plot by event type and comparison
    p3b <- ggplot(sig_events, aes(x = event_type, y = IncLevelDifference, fill = event_type)) +
        geom_boxplot(alpha = 0.7, outlier.shape = 21) +
        facet_wrap(~ comparison, ncol = 2) +
        scale_fill_manual(values = event_colors) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            title = "deltaPSI Distribution by Event Type",
            x = "Event Type",
            y = "deltaPSI",
            fill = "Event Type"
        )

    ggsave(file.path(OUTPUT_DIR, "dpsi_boxplot_by_type.pdf"), p3b, width = 12, height = 10)

    # 3c. Histogram with directionality
    p3c <- ggplot(sig_events, aes(x = IncLevelDifference, fill = direction)) +
        geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
        facet_wrap(~ comparison, ncol = 2, scales = "free_y") +
        scale_fill_manual(values = c("Inclusion" = "#e74c3c", "Exclusion" = "#3498db")) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
        theme_bw() +
        labs(
            title = "deltaPSI Distribution: Inclusion vs Exclusion",
            x = "deltaPSI",
            y = "Count",
            fill = "Direction"
        )

    ggsave(file.path(OUTPUT_DIR, "dpsi_histogram_direction.pdf"), p3c, width = 12, height = 10)
}

#===============================================================================
# 4. Directionality Analysis
#===============================================================================

cat("Creating directionality plots...\n")

if (nrow(sig_events) > 0) {
    direction_summary <- sig_events %>%
        group_by(comparison, event_type, direction) %>%
        summarise(count = n(), .groups = "drop")

    # Directionality bar plot
    p4 <- ggplot(direction_summary, aes(x = event_type, y = count, fill = direction)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~ comparison, ncol = 2) +
        scale_fill_manual(values = c("Inclusion" = "#e74c3c", "Exclusion" = "#3498db")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            title = "Splicing Directionality: Inclusion vs Exclusion",
            subtitle = "Positive deltaPSI = Inclusion, Negative deltaPSI = Exclusion",
            x = "Event Type",
            y = "Number of Events",
            fill = "Direction"
        )

    ggsave(file.path(OUTPUT_DIR, "directionality_barplot.pdf"), p4, width = 12, height = 10)
    ggsave(file.path(OUTPUT_DIR, "directionality_barplot.png"), p4, width = 12, height = 10, dpi = 300)

    # Save directionality summary
    direction_wide <- direction_summary %>%
        pivot_wider(names_from = direction, values_from = count, values_fill = 0) %>%
        mutate(
            total = Inclusion + Exclusion,
            inclusion_ratio = round(Inclusion / total, 3)
        )

    write.csv(direction_wide, file.path(OUTPUT_DIR, "directionality_summary.csv"), row.names = FALSE)
}

#===============================================================================
# 5. Venn Diagram of Overlapping Genes
#===============================================================================

cat("Creating Venn diagram...\n")

# Get significant genes for each comparison
sig_genes_list <- list()
for (comp in comparisons) {
    genes <- combined_df %>%
        filter(comparison == comp, significant) %>%
        pull(geneSymbol) %>%
        unique()
    sig_genes_list[[comp]] <- genes
    cat("  ", comp, ":", length(genes), "genes with significant splicing\n")
}

# Create Venn diagram (supports 2-4 sets)
if (length(sig_genes_list) >= 2 && all(sapply(sig_genes_list, length) > 0)) {
    # Use shortened names for display
    short_names <- gsub("_vs_", "\nvs\n", names(sig_genes_list))
    names(sig_genes_list) <- short_names

    # Generate Venn diagram
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

    pdf(file.path(OUTPUT_DIR, "venn_diagram_genes.pdf"), width = 12, height = 10)

    if (length(sig_genes_list) == 2) {
        venn_plot <- draw.pairwise.venn(
            area1 = length(sig_genes_list[[1]]),
            area2 = length(sig_genes_list[[2]]),
            cross.area = length(intersect(sig_genes_list[[1]], sig_genes_list[[2]])),
            category = names(sig_genes_list),
            fill = c("#E41A1C", "#377EB8"),
            alpha = 0.5,
            cex = 1.5,
            cat.cex = 1.2
        )
    } else if (length(sig_genes_list) == 3) {
        venn_plot <- draw.triple.venn(
            area1 = length(sig_genes_list[[1]]),
            area2 = length(sig_genes_list[[2]]),
            area3 = length(sig_genes_list[[3]]),
            n12 = length(intersect(sig_genes_list[[1]], sig_genes_list[[2]])),
            n23 = length(intersect(sig_genes_list[[2]], sig_genes_list[[3]])),
            n13 = length(intersect(sig_genes_list[[1]], sig_genes_list[[3]])),
            n123 = length(Reduce(intersect, sig_genes_list)),
            category = names(sig_genes_list),
            fill = c("#E41A1C", "#377EB8", "#4DAF4A"),
            alpha = 0.5,
            cex = 1.5,
            cat.cex = 1.0
        )
    } else if (length(sig_genes_list) == 4) {
        venn_plot <- draw.quad.venn(
            area1 = length(sig_genes_list[[1]]),
            area2 = length(sig_genes_list[[2]]),
            area3 = length(sig_genes_list[[3]]),
            area4 = length(sig_genes_list[[4]]),
            n12 = length(intersect(sig_genes_list[[1]], sig_genes_list[[2]])),
            n13 = length(intersect(sig_genes_list[[1]], sig_genes_list[[3]])),
            n14 = length(intersect(sig_genes_list[[1]], sig_genes_list[[4]])),
            n23 = length(intersect(sig_genes_list[[2]], sig_genes_list[[3]])),
            n24 = length(intersect(sig_genes_list[[2]], sig_genes_list[[4]])),
            n34 = length(intersect(sig_genes_list[[3]], sig_genes_list[[4]])),
            n123 = length(Reduce(intersect, sig_genes_list[1:3])),
            n124 = length(Reduce(intersect, sig_genes_list[c(1,2,4)])),
            n134 = length(Reduce(intersect, sig_genes_list[c(1,3,4)])),
            n234 = length(Reduce(intersect, sig_genes_list[2:4])),
            n1234 = length(Reduce(intersect, sig_genes_list)),
            category = names(sig_genes_list),
            fill = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
            alpha = 0.5,
            cex = 1.2,
            cat.cex = 0.9
        )
    }

    grid.draw(venn_plot)
    grid.text("Genes with Significant Splicing Events\n(Overlap across comparisons)",
              x = 0.5, y = 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
    dev.off()
}

#===============================================================================
# 6. Summary Tables
#===============================================================================

cat("Creating summary tables...\n")

# 6a. Event type summary
event_summary <- summary_df %>%
    group_by(event_type) %>%
    summarise(
        total_all_comparisons = sum(total),
        significant_all_comparisons = sum(significant),
        pct_significant = round(100 * significant_all_comparisons / total_all_comparisons, 1),
        .groups = "drop"
    ) %>%
    arrange(desc(significant_all_comparisons))

write.csv(event_summary, file.path(OUTPUT_DIR, "event_type_summary.csv"), row.names = FALSE)

# 6b. Comparison summary
comparison_summary <- combined_df %>%
    group_by(comparison) %>%
    summarise(
        total_events = n(),
        significant_events = sum(significant),
        affected_genes = n_distinct(geneSymbol[significant]),
        mean_abs_dpsi = round(mean(abs(IncLevelDifference[significant])), 3),
        .groups = "drop"
    )

write.csv(comparison_summary, file.path(OUTPUT_DIR, "comparison_summary.csv"), row.names = FALSE)

# 6c. Top significant events with IGV coordinates
top_events <- sig_events %>%
    arrange(desc(abs(IncLevelDifference))) %>%
    mutate(
        IGV_coord = paste0(chr, ":",
                          pmax(0, as.integer(exonStart_0base) - 500), "-",
                          as.integer(exonEnd) + 500)
    ) %>%
    select(comparison, event_type, geneSymbol, GeneID, chr, strand,
           IncLevelDifference, min_total_reads, direction, IGV_coord) %>%
    head(200)

write.csv(top_events, file.path(OUTPUT_DIR, "top_significant_events_with_IGV.csv"), row.names = FALSE)

# 6d. Gene-level summary
gene_summary <- sig_events %>%
    group_by(geneSymbol) %>%
    summarise(
        n_events = n(),
        comparisons = paste(unique(comparison), collapse = "; "),
        event_types = paste(unique(event_type), collapse = "; "),
        max_abs_dpsi = round(max(abs(IncLevelDifference)), 3),
        mean_dpsi = round(mean(IncLevelDifference), 3),
        .groups = "drop"
    ) %>%
    arrange(desc(n_events))

write.csv(gene_summary, file.path(OUTPUT_DIR, "gene_level_summary.csv"), row.names = FALSE)

#===============================================================================
# 7. Overlap Heatmap Between Comparisons
#===============================================================================

cat("Creating overlap heatmap...\n")

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

# Save overlap matrix
write.csv(as.data.frame(overlap_matrix), file.path(OUTPUT_DIR, "gene_overlap_matrix.csv"))

# Create heatmap
pdf(file.path(OUTPUT_DIR, "comparison_overlap_heatmap.pdf"), width = 10, height = 8)
pheatmap(overlap_matrix,
         display_numbers = TRUE,
         number_format = "%d",
         main = "Overlap of Genes with Significant Splicing Events",
         color = colorRampPalette(brewer.pal(9, "Blues"))(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize = 10,
         fontsize_number = 12)
dev.off()

#===============================================================================
# Summary Report
#===============================================================================

cat("\n=== Summary Report ===\n")
cat("Filtering criteria applied:\n")
cat("  - |deltaPSI| >", PSI_THRESHOLD, "\n")
cat("  - Minimum junction reads >=", MIN_READS, "\n")
cat("\n")
print(comparison_summary)

cat("\n=== IMPORTANT CAVEATS ===\n")
cat("This analysis was performed on MONOPLICATE data (n=1 per condition).\n")
cat("Results should be interpreted as EXPLORATORY only.\n")
cat("False positive rate is expected to be HIGH.\n")
cat("Key findings should be validated experimentally (RT-PCR, etc.)\n")

cat("\n=== Visualization Complete ===\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Files generated:\n")
cat("  - event_distribution_*.pdf/png  : Event type distributions\n")
cat("  - pie_chart_*.pdf               : Pie charts per comparison\n")
cat("  - dpsi_*.pdf/png                : deltaPSI distribution plots\n")
cat("  - directionality_*.pdf/csv      : Inclusion vs exclusion analysis\n")
cat("  - venn_diagram_genes.pdf        : Gene overlap Venn diagram\n")
cat("  - *_summary.csv                 : Summary tables\n")
cat("  - top_significant_events_with_IGV.csv : Top events with IGV coordinates\n")
cat("\nTimestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
