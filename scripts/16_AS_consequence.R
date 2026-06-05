#!/usr/bin/env Rscript
# =============================================================================
# 16_AS_consequence.R  --  Functional consequence of each splicing event
# =============================================================================
#
# A splicing event only MEANS something through what it does to the transcript.
# This step annotates every filtered event with its likely protein consequence,
# so "gene X is spliced" becomes "exon skipping in X is frame-shifting, overlaps
# the CDS, and is more included in the mutant -> likely PTC -> NMD".
#
# THREE ANNOTATIONS (coordinate-based, no fragile downloads required):
#   1. FRAME IMPACT      length of the affected segment %% 3
#        SE  : cassette exon length            (frame-shifting if not 3n)
#        RI  : retained intron length          (almost always frame-shifting)
#        A5SS/A3SS: |long - short| exon length (frame-shifting if not 3n)
#   2. CDS vs UTR        does the affected segment overlap any CDS in the GTF?
#        no overlap  -> UTR / non-coding: affects stability/translation, not
#                       protein sequence (e.g. 3'UTR -> miRNA/decay).
#        overlap     -> coding: frame impact + domain consequences apply.
#   3. NMD PREDICTION    coding & frame-shifting & MORE-included-in-mutant
#        => introduces a premature stop in the mutant => NMD => expected
#           DOWN-regulation. We test this against the DEG log2FC (a built-in
#           sanity check of the NMD interpretation).
#
# Optional: if the `maser` package is installed, per-transcript plots for the
# top concordant genes are emitted too (purely illustrative).
#
# n = 1 CAVEAT: consequences are predicted from annotation, not measured.
# Treat as hypotheses; the DEG correlation is corroborating, not proof.
#
# INPUTS:
#   results/05_splicing/<comp>/{SE,RI,A5SS,A3SS}.MATS.JC.txt   (raw rMATS)
#   results/04_DEG/DEG_<comp>.csv                              (for NMD check)
#   GTF (GRCm39) for CDS features
# OUTPUTS -> results/09_AS_interpretation/consequence/
#   <comp>_event_consequences.csv
#   <comp>_consequence_summary.txt
#   <comp>_nmd_vs_deg.png
#   consequence_overview.txt
#
# USAGE (cluster):
#   Rscript 16_AS_consequence.R \
#     --base-dir /path/to/90-1265107649 \
#     --gtf /beegfs/.../refdata-gex-GRCm39-2024-A/genes/genes.gtf
# =============================================================================

suppressWarnings(suppressMessages({
  library(rtracklayer)
  library(GenomicRanges)
}))

# Minimal base-R arg parsing (avoids an optparse dependency on the cluster).
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- which(args == flag)
  if (length(i) == 1L && i < length(args)) return(args[i + 1L])
  default
}

base_dir  <- get_arg("--base-dir",
                     "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649")
gtf_path  <- get_arg("--gtf",
                     "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf")
min_dpsi  <- as.numeric(get_arg("--min-dpsi", "0.10"))
min_reads <- as.integer(get_arg("--min-reads", "20"))

comparisons <- c("EMX1_wt_vs_mut", "Nestin_wt_vs_mut")
splic_dir   <- file.path(base_dir, "results", "05_splicing")
deg_dir     <- file.path(base_dir, "results", "04_DEG")
out_dir     <- file.path(base_dir, "results", "09_AS_interpretation", "consequence")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- CDS GRanges from the GTF ------------------------------------------------
message("[16] importing CDS features from GTF ...")
gtf <- rtracklayer::import(gtf_path)
cds <- gtf[gtf$type == "CDS"]
cds <- reduce(cds)               # collapse to genomic CDS footprint
message(sprintf("[16] %d reduced CDS ranges.", length(cds)))

first_num <- function(x) {
  # IncLevel can be comma-joined across replicates; n=1 -> single value
  suppressWarnings(as.numeric(sapply(strsplit(as.character(x), ","), `[`, 1)))
}

affected_region <- function(df, event) {
  # returns data.frame(start, end, seg_len) on the + strand genomic coords
  if (event == "SE") {
    data.frame(start = df$exonStart_0base, end = df$exonEnd,
               seg_len = df$exonEnd - df$exonStart_0base)
  } else if (event == "RI") {
    data.frame(start = df$upstreamEE, end = df$downstreamES,
               seg_len = df$downstreamES - df$upstreamEE)
  } else { # A5SS / A3SS
    long_len  <- df$longExonEnd - df$longExonStart_0base
    short_len <- df$shortEE - df$shortES
    data.frame(start = pmin(df$shortES, df$longExonStart_0base),
               end   = pmax(df$shortEE, df$longExonEnd),
               seg_len = abs(long_len - short_len))
  }
}

overview <- c("AS interpretation -- step 16 (functional consequence)",
              paste(rep("=", 60), collapse = ""),
              "Frame impact + CDS/UTR + NMD prediction; NMD checked vs DEG.",
              "Predicted from annotation (n=1) -- hypotheses, not measurements.", "")

for (comp in comparisons) {
  message(sprintf("[16] %s ...", comp))
  events_all <- list()
  for (event in c("SE", "RI", "A5SS", "A3SS")) {
    fp <- file.path(splic_dir, comp, paste0(event, ".MATS.JC.txt"))
    if (!file.exists(fp)) next
    df <- read.delim(fp, stringsAsFactors = FALSE)

    r1 <- df$IJC_SAMPLE_1 + df$SJC_SAMPLE_1
    r2 <- df$IJC_SAMPLE_2 + df$SJC_SAMPLE_2
    ild <- suppressWarnings(as.numeric(df$IncLevelDifference))
    keep <- which(abs(ild) >= min_dpsi & r1 >= min_reads & r2 >= min_reads &
                    !is.na(ild))
    if (length(keep) == 0) next
    df <- df[keep, , drop = FALSE]
    ild <- ild[keep]                          # subset to match the kept rows

    reg <- affected_region(df, event)
    gr <- GRanges(seqnames = df$chr,
                  ranges = IRanges(start = reg$start + 1L, end = reg$end),
                  strand = df$strand)
    in_cds <- overlapsAny(gr, cds, ignore.strand = TRUE)

    frame_shift <- (reg$seg_len %% 3L) != 0L
    dpsi_mut <- -ild                          # >0 => more inclusion in mutant
    direction <- ifelse(dpsi_mut > 0, "inclusion_in_mut", "exclusion_in_mut")

    # NMD logic: a frame-shifting CDS segment introduces a PTC when included.
    # If it is MORE included in the mutant -> PTC gained in mutant -> NMD ->
    # expect lower expression. If MORE excluded in mutant -> PTC removed.
    nmd_pred <- ifelse(
      in_cds & frame_shift & dpsi_mut > 0, "PTC_gained_in_mut(->down)",
      ifelse(in_cds & frame_shift & dpsi_mut < 0, "PTC_lost_in_mut(->up)",
             ifelse(in_cds & !frame_shift, "frame_preserving_coding",
                    "UTR_or_noncoding")))

    events_all[[event]] <- data.frame(
      comparison = comp, event_type = event,
      GeneID = df$GeneID, geneSymbol = df$geneSymbol,
      chr = df$chr, strand = df$strand,
      seg_len = reg$seg_len, frame_shifting = frame_shift,
      in_cds = in_cds, region_class = ifelse(in_cds, "CDS", "UTR/noncoding"),
      IncLevelDifference = ild, dPSI_mut = dpsi_mut, direction = direction,
      nmd_prediction = nmd_pred,
      stringsAsFactors = FALSE)
  }
  if (length(events_all) == 0) { overview <- c(overview, paste0(comp, ": none")); next }
  ev <- do.call(rbind, events_all)
  ev$geneSymbol <- as.character(ev$geneSymbol)

  # ---- NMD vs DEG sanity check ----------------------------------------------
  deg_fp <- file.path(deg_dir, paste0("DEG_", comp, ".csv"))
  nmd_line <- "DEG file not found -- NMD-vs-expression check skipped."
  if (file.exists(deg_fp)) {
    deg <- read.csv(deg_fp, stringsAsFactors = FALSE)
    # log2FoldChange is mutant-relative (sample2/sample1); <0 = down in mutant
    ev$gene_log2FC <- deg$log2FoldChange[match(ev$GeneID, deg$gene_id)]
    ptc_up   <- ev$gene_log2FC[ev$nmd_prediction == "PTC_gained_in_mut(->down)"]
    others   <- ev$gene_log2FC[ev$nmd_prediction == "frame_preserving_coding"]
    ptc_up <- ptc_up[is.finite(ptc_up)]; others <- others[is.finite(others)]
    if (length(ptc_up) >= 5 && length(others) >= 5) {
      wt <- suppressWarnings(wilcox.test(ptc_up, others, alternative = "less"))
      nmd_line <- sprintf(
        "NMD check: PTC-gained events median log2FC=%.2f (n=%d) vs frame-preserving median=%.2f (n=%d); Wilcoxon(less) p=%.3g %s",
        median(ptc_up), length(ptc_up), median(others), length(others),
        wt$p.value,
        ifelse(wt$p.value < 0.05, "-> consistent with NMD.",
               "-> not significant (n=1; expected to be weak)."))
      png(file.path(out_dir, paste0(comp, "_nmd_vs_deg.png")),
          width = 1100, height = 800, res = 150)
      boxplot(list(`PTC gained\n(->NMD)` = ptc_up,
                   `frame-preserving\ncoding` = others),
              ylab = "gene log2FC (mutant vs wt)",
              main = paste0("NMD prediction vs expression: ", comp,
                            "\n[exploratory, n=1]"),
              col = c("salmon", "lightblue"))
      abline(h = 0, lty = 2)
      dev.off()
    }
  }

  write.csv(ev, file.path(out_dir, paste0(comp, "_event_consequences.csv")),
            row.names = FALSE)

  tab <- table(ev$region_class, ev$frame_shifting)
  summ <- c(
    sprintf("%s: %d filtered events", comp, nrow(ev)),
    sprintf("  CDS events: %d (%.0f%%); UTR/noncoding: %d",
            sum(ev$in_cds), 100 * mean(ev$in_cds), sum(!ev$in_cds)),
    sprintf("  frame-shifting (all): %d (%.0f%%)",
            sum(ev$frame_shifting), 100 * mean(ev$frame_shifting)),
    sprintf("  predicted PTC-gained-in-mutant (NMD candidates): %d",
            sum(ev$nmd_prediction == "PTC_gained_in_mut(->down)")),
    paste0("  ", nmd_line))
  writeLines(c(paste0(comp, " consequence summary"), summ),
             file.path(out_dir, paste0(comp, "_consequence_summary.txt")))
  overview <- c(overview, summ, "")
}

writeLines(overview, file.path(out_dir, "consequence_overview.txt"))
cat(paste(overview, collapse = "\n"), "\n")
cat(sprintf("\n[16] outputs in %s\n", out_dir))
