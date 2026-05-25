#!/usr/bin/env Rscript

#===============================================================================
# SCRIPT: 5_DEG_analysis.R
# PURPOSE: Differential Gene Expression Analysis for RiboTag Data
#
# DESCRIPTION:
# Performs fold-change based differential expression analysis using DESeq2
# for normalization. Since samples are monoplicates (n=1), we rely on
# fold-change cutoffs rather than statistical p-values.
#
# COMPARISONS:
# 1. EMX1: wt vs mut
# 2. Nestin: wt vs mut
# 3. WT: EMX1 vs Nestin
# 4. MUT: EMX1 vs Nestin
#
# IMPORTANT - MONOPLICATE DESIGN LIMITATIONS:
# This experiment has n=1 per condition (no biological replicates).
# Statistical significance testing (p-values) is NOT possible.
# Results are based on fold-change cutoffs only (|log2FC| > 1).
# Interpretation should be exploratory, not definitive.
# Validation of key findings is strongly recommended.
#
# LIBRARY STRANDEDNESS:
# Verified as UNSTRANDED (Forward/Reverse ratio ≈ 0.97)
# Using column 2 (unstranded) from STAR ReadsPerGene.out.tab
#
# USAGE:
# Rscript 5_DEG_analysis.R
#===============================================================================

# Load libraries
suppressPackageStartupMessages({
    library(DESeq2)
    library(tidyverse)
    library(pheatmap)
    library(RColorBrewer)
    library(ggrepel)
})

cat("=== Starting DEG Analysis ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Set paths
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
ALIGNED_DIR <- file.path(BASE_DIR, "results/02_aligned")
OUTPUT_DIR <- file.path(BASE_DIR, "results/04_DEG")

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# 1. Load count data from STAR output
#-------------------------------------------------------------------------------
cat("=== Step 1: Loading count data ===\n")

# Sample metadata
sample_info <- data.frame(
    sample_id = c("1", "2", "3", "4"),
    sample_name = c("EMX1_hippo_wt", "EMX1_hippo_mut", "Nestin_hippo_wt", "Nestin_hippo_mut"),
    cre_driver = c("EMX1", "EMX1", "Nestin", "Nestin"),
    genotype = c("wt", "mut", "wt", "mut"),
    stringsAsFactors = FALSE
)

# Read count files from STAR
count_files <- file.path(ALIGNED_DIR, sample_info$sample_id,
                         paste0(sample_info$sample_id, "_ReadsPerGene.out.tab"))

# Check files exist
for (f in count_files) {
    if (!file.exists(f)) {
        stop("Count file not found: ", f)
    }
}

# Read and merge count data
# STAR ReadsPerGene.out.tab has 4 columns: gene_id, unstranded, forward, reverse
# We'll use column 2 (unstranded) - adjust if your library is stranded
read_star_counts <- function(file, sample_name) {
    df <- read.table(file, header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("gene_id", "unstranded", "forward", "reverse"))
    # Skip first 4 rows (N_unmapped, N_multimapping, N_noFeature, N_ambiguous)
    df <- df[!grepl("^N_", df$gene_id), ]
    df <- df[, c("gene_id", "unstranded")]
    colnames(df)[2] <- sample_name
    return(df)
}

# Merge all count files
counts_list <- lapply(seq_len(nrow(sample_info)), function(i) {
    read_star_counts(count_files[i], sample_info$sample_name[i])
})

counts_df <- counts_list[[1]]
for (i in 2:length(counts_list)) {
    counts_df <- merge(counts_df, counts_list[[i]], by = "gene_id")
}

# Set gene_id as rownames
rownames(counts_df) <- counts_df$gene_id
counts_df$gene_id <- NULL
counts_matrix <- as.matrix(counts_df)

cat("Loaded counts for", nrow(counts_matrix), "genes across", ncol(counts_matrix), "samples\n")

# Save raw counts
write.csv(counts_df, file.path(OUTPUT_DIR, "raw_counts.csv"))

#-------------------------------------------------------------------------------
# 2. DESeq2 Normalization
#-------------------------------------------------------------------------------
cat("\n=== Step 2: DESeq2 Normalization ===\n")

# Create DESeq2 dataset (minimal design for normalization only)
rownames(sample_info) <- sample_info$sample_name
dds <- DESeqDataSetFromMatrix(
    countData = counts_matrix,
    colData = sample_info,
    design = ~ 1  # Minimal design for normalization
)

# Filter low-count genes (at least 10 counts in at least 2 samples)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

# Estimate size factors for normalization
dds <- estimateSizeFactors(dds)
cat("Size factors:", sizeFactors(dds), "\n")

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file.path(OUTPUT_DIR, "normalized_counts.csv"))

# VST transformation for visualization
vst_counts <- vst(dds, blind = TRUE)

#-------------------------------------------------------------------------------
# 3. Sample QC Plots
#-------------------------------------------------------------------------------
cat("\n=== Step 3: Generating QC plots ===\n")

# PCA plot
pca_data <- plotPCA(vst_counts, intgroup = c("cre_driver", "genotype"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = cre_driver, shape = genotype)) +
    geom_point(size = 5) +
    geom_text_repel(aes(label = name), size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_bw() +
    ggtitle("PCA of RiboTag Samples") +
    scale_color_brewer(palette = "Set1")

ggsave(file.path(OUTPUT_DIR, "PCA_plot.pdf"), pca_plot, width = 8, height = 6)
ggsave(file.path(OUTPUT_DIR, "PCA_plot.png"), pca_plot, width = 8, height = 6, dpi = 300)

# Sample correlation heatmap
cor_matrix <- cor(normalized_counts, method = "spearman")
pdf(file.path(OUTPUT_DIR, "sample_correlation_heatmap.pdf"), width = 8, height = 7)
pheatmap(cor_matrix,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         display_numbers = TRUE,
         number_format = "%.3f",
         main = "Sample Correlation (Spearman)")
dev.off()

#-------------------------------------------------------------------------------
# 4. Fold-Change Analysis Function
#-------------------------------------------------------------------------------
cat("\n=== Step 4: Computing fold changes ===\n")

compute_fold_change <- function(norm_counts, sample1, sample2, comparison_name,
                                fc_cutoff = 1, output_dir = OUTPUT_DIR) {

    cat("\nProcessing comparison:", comparison_name, "\n")
    cat("  ", sample1, "vs", sample2, "\n")

    # Get counts for the two samples
    counts1 <- norm_counts[, sample1]
    counts2 <- norm_counts[, sample2]

    # Calculate log2 fold change (add pseudocount to avoid log(0))
    pseudocount <- 1
    log2FC <- log2((counts2 + pseudocount) / (counts1 + pseudocount))

    # Calculate mean expression
    baseMean <- (counts1 + counts2) / 2

    # Create results dataframe
    results <- data.frame(
        gene_id = rownames(norm_counts),
        baseMean = baseMean,
        counts_sample1 = counts1,
        counts_sample2 = counts2,
        log2FoldChange = log2FC,
        absLog2FC = abs(log2FC),
        stringsAsFactors = FALSE
    )

    # Classify genes
    results$regulation <- "NS"  # Not significant
    results$regulation[results$log2FoldChange > fc_cutoff & results$baseMean > 10] <- "UP"
    results$regulation[results$log2FoldChange < -fc_cutoff & results$baseMean > 10] <- "DOWN"

    # Sort by absolute fold change
    results <- results[order(-results$absLog2FC), ]

    # Summary
    cat("  Upregulated (log2FC >", fc_cutoff, "):", sum(results$regulation == "UP"), "\n")
    cat("  Downregulated (log2FC <", -fc_cutoff, "):", sum(results$regulation == "DOWN"), "\n")

    # Save results
    output_file <- file.path(output_dir, paste0("DEG_", comparison_name, ".csv"))
    write.csv(results, output_file, row.names = FALSE)

    # Save significant genes only
    sig_results <- results[results$regulation != "NS", ]
    sig_file <- file.path(output_dir, paste0("DEG_", comparison_name, "_significant.csv"))
    write.csv(sig_results, sig_file, row.names = FALSE)

    # MA plot
    ma_plot <- ggplot(results, aes(x = log10(baseMean + 1), y = log2FoldChange, color = regulation)) +
        geom_point(alpha = 0.5, size = 1) +
        scale_color_manual(values = c("DOWN" = "blue", "NS" = "gray", "UP" = "red")) +
        geom_hline(yintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black") +
        geom_hline(yintercept = 0, color = "black") +
        theme_bw() +
        labs(title = paste0("MA Plot: ", comparison_name),
             subtitle = paste0("UP: ", sum(results$regulation == "UP"),
                             " | DOWN: ", sum(results$regulation == "DOWN")),
             x = "log10(Mean Expression + 1)",
             y = "log2 Fold Change") +
        theme(legend.position = "bottom")

    ggsave(file.path(output_dir, paste0("MA_plot_", comparison_name, ".pdf")),
           ma_plot, width = 8, height = 6)
    ggsave(file.path(output_dir, paste0("MA_plot_", comparison_name, ".png")),
           ma_plot, width = 8, height = 6, dpi = 300)

    return(results)
}

#-------------------------------------------------------------------------------
# 5. Run All Comparisons
#-------------------------------------------------------------------------------
cat("\n=== Step 5: Running all comparisons ===\n")

# Comparison 1: EMX1 wt vs mut
deg_emx1 <- compute_fold_change(normalized_counts,
                                 "EMX1_hippo_wt", "EMX1_hippo_mut",
                                 "EMX1_wt_vs_mut")

# Comparison 2: Nestin wt vs mut
deg_nestin <- compute_fold_change(normalized_counts,
                                   "Nestin_hippo_wt", "Nestin_hippo_mut",
                                   "Nestin_wt_vs_mut")

# Comparison 3: WT - EMX1 vs Nestin
deg_wt_cre <- compute_fold_change(normalized_counts,
                                   "EMX1_hippo_wt", "Nestin_hippo_wt",
                                   "WT_EMX1_vs_Nestin")

# Comparison 4: MUT - EMX1 vs Nestin
deg_mut_cre <- compute_fold_change(normalized_counts,
                                    "EMX1_hippo_mut", "Nestin_hippo_mut",
                                    "MUT_EMX1_vs_Nestin")

#-------------------------------------------------------------------------------
# 6. Combined Visualization
#-------------------------------------------------------------------------------
cat("\n=== Step 6: Creating combined visualizations ===\n")

# Heatmap of top variable genes
top_var_genes <- head(order(rowVars(as.matrix(normalized_counts)), decreasing = TRUE), 50)
top_genes_matrix <- normalized_counts[top_var_genes, ]

# Z-score normalization for heatmap
top_genes_zscore <- t(scale(t(top_genes_matrix)))

# Annotation
annotation_col <- data.frame(
    Cre_Driver = sample_info$cre_driver,
    Genotype = sample_info$genotype,
    row.names = sample_info$sample_name
)

pdf(file.path(OUTPUT_DIR, "top50_variable_genes_heatmap.pdf"), width = 10, height = 12)
pheatmap(top_genes_zscore,
         annotation_col = annotation_col,
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         fontsize_row = 6,
         main = "Top 50 Most Variable Genes (Z-score)")
dev.off()

#-------------------------------------------------------------------------------
# 6b. Heatmap with Gene Names (instead of Gene IDs)
#-------------------------------------------------------------------------------
cat("\n=== Step 6b: Creating heatmap with gene names ===\n")

# Load GTF to get gene ID to gene name mapping (Mouse GRCm39)
GTF_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"

# Read GTF and extract gene_id to gene_name mapping
cat("Loading gene annotations from GTF...\n")
gtf_genes <- tryCatch({
    # Read only gene lines from GTF
    gtf_lines <- readLines(GTF_FILE)
    gene_lines <- gtf_lines[grepl("\tgene\t", gtf_lines)]

    # Extract gene_id and gene_name
    gene_id <- gsub('.*gene_id "([^"]+)".*', '\\1', gene_lines)
    gene_name <- gsub('.*gene_name "([^"]+)".*', '\\1', gene_lines)

    # Create mapping dataframe
    data.frame(gene_id = gene_id, gene_name = gene_name, stringsAsFactors = FALSE)
}, error = function(e) {
    cat("Warning: Could not load GTF file. Using gene IDs as names.\n")
    NULL
})

if (!is.null(gtf_genes)) {
    cat("Loaded", nrow(gtf_genes), "gene annotations\n")

    # Create gene ID to name mapping (handle version numbers in ENSEMBL IDs)
    # Gene IDs in counts may have version (ENSG00000000003.15) or not (ENSG00000000003)
    gtf_genes$gene_id_no_version <- gsub("\\..*", "", gtf_genes$gene_id)

    # Map gene names to our top genes
    top_gene_ids <- rownames(top_genes_zscore)
    top_gene_ids_no_version <- gsub("\\..*", "", top_gene_ids)

    # Create mapping
    gene_name_map <- gtf_genes$gene_name[match(top_gene_ids_no_version, gtf_genes$gene_id_no_version)]

    # Handle unmapped genes (keep original ID if no name found)
    gene_name_map[is.na(gene_name_map)] <- top_gene_ids[is.na(gene_name_map)]

    # Handle duplicates by appending a suffix
    if (any(duplicated(gene_name_map))) {
        dup_names <- gene_name_map[duplicated(gene_name_map)]
        for (dup in unique(dup_names)) {
            idx <- which(gene_name_map == dup)
            gene_name_map[idx] <- paste0(dup, "_", seq_along(idx))
        }
    }

    # Create heatmap matrix with gene names
    top_genes_zscore_named <- top_genes_zscore
    rownames(top_genes_zscore_named) <- gene_name_map

    # Generate heatmap with gene names
    pdf(file.path(OUTPUT_DIR, "top50_variable_genes_heatmap_genenames.pdf"), width = 10, height = 12)
    pheatmap(top_genes_zscore_named,
             annotation_col = annotation_col,
             clustering_method = "ward.D2",
             show_rownames = TRUE,
             fontsize_row = 6,
             main = "Top 50 Most Variable Genes (Z-score) - Gene Names")
    dev.off()

    # Also save PNG version
    png(file.path(OUTPUT_DIR, "top50_variable_genes_heatmap_genenames.png"),
        width = 10, height = 12, units = "in", res = 300)
    pheatmap(top_genes_zscore_named,
             annotation_col = annotation_col,
             clustering_method = "ward.D2",
             show_rownames = TRUE,
             fontsize_row = 6,
             main = "Top 50 Most Variable Genes (Z-score) - Gene Names")
    dev.off()

    cat("Created heatmap with gene names: top50_variable_genes_heatmap_genenames.pdf\n")

    # Save gene ID to name mapping for the top 50 genes
    top_genes_mapping <- data.frame(
        gene_id = top_gene_ids,
        gene_name = gene_name_map,
        stringsAsFactors = FALSE
    )
    write.csv(top_genes_mapping, file.path(OUTPUT_DIR, "top50_genes_id_to_name.csv"), row.names = FALSE)
}

# Venn diagram data - genes affected in both EMX1 and Nestin
emx1_up <- deg_emx1$gene_id[deg_emx1$regulation == "UP"]
emx1_down <- deg_emx1$gene_id[deg_emx1$regulation == "DOWN"]
nestin_up <- deg_nestin$gene_id[deg_nestin$regulation == "UP"]
nestin_down <- deg_nestin$gene_id[deg_nestin$regulation == "DOWN"]

# Common genes
common_up <- intersect(emx1_up, nestin_up)
common_down <- intersect(emx1_down, nestin_down)

cat("\nGenes upregulated in both EMX1 and Nestin (mut vs wt):", length(common_up), "\n")
cat("Genes downregulated in both EMX1 and Nestin (mut vs wt):", length(common_down), "\n")

# Save common genes
if (length(common_up) > 0) {
    write.csv(data.frame(gene_id = common_up),
              file.path(OUTPUT_DIR, "common_upregulated_genes.csv"), row.names = FALSE)
}
if (length(common_down) > 0) {
    write.csv(data.frame(gene_id = common_down),
              file.path(OUTPUT_DIR, "common_downregulated_genes.csv"), row.names = FALSE)
}

#-------------------------------------------------------------------------------
# 7. Summary Report
#-------------------------------------------------------------------------------
cat("\n=== Step 7: Generating summary report ===\n")

summary_df <- data.frame(
    Comparison = c("EMX1_wt_vs_mut", "Nestin_wt_vs_mut", "WT_EMX1_vs_Nestin", "MUT_EMX1_vs_Nestin"),
    Upregulated = c(sum(deg_emx1$regulation == "UP"),
                    sum(deg_nestin$regulation == "UP"),
                    sum(deg_wt_cre$regulation == "UP"),
                    sum(deg_mut_cre$regulation == "UP")),
    Downregulated = c(sum(deg_emx1$regulation == "DOWN"),
                      sum(deg_nestin$regulation == "DOWN"),
                      sum(deg_wt_cre$regulation == "DOWN"),
                      sum(deg_mut_cre$regulation == "DOWN"))
)

write.csv(summary_df, file.path(OUTPUT_DIR, "DEG_summary.csv"), row.names = FALSE)

cat("\n=== DEG Analysis Summary ===\n")
print(summary_df)

cat("\n=== DEG Analysis Complete ===\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# Save session info
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))
