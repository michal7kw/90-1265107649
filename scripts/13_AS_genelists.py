#!/usr/bin/env python3
"""
================================================================================
13_AS_genelists.py  --  Foundation step for AS interpretation (steps 13-17)
================================================================================

Turns the raw rMATS differential-splicing output into the curated gene lists,
the cross-lineage concordant set, and the background universe that every
downstream interpretation step (14 enrichment, 15 rMAPS, 16 consequence,
17 integration) consumes.

SCOPE (per the approved plan):
  - Only the two MUTATION-effect comparisons: EMX1_wt_vs_mut and Nestin_wt_vs_mut.
  - Event types SE, A5SS, A3SS, RI (MXE dropped -- 0 significant in this dataset).

WHY A CUSTOM FILTER (n = 1 per condition):
  rMATS PValue/FDR are NOT trustworthy with monoplicate data. We therefore gate
  on effect size + read support, exactly mirroring the filter the upstream
  pipeline (6_differential_splicing.sh) already used for its summary:
        |IncLevelDifference| >= MIN_DPSI   (default 0.10)
        AND  IJC+SJC >= MIN_READS in BOTH samples (default 20)
  PValue/FDR are carried through as columns but never used as the gate.

DIRECTION CONVENTION:
  rMATS reports IncLevelDifference = IncLevel1 - IncLevel2 = (WT) - (MUT) for
  these comparisons (b1 = wt sample, b2 = mut sample in 6_differential_splicing.sh).
  We therefore define
        dPSI_mut = IncLevel2 - IncLevel1 = -IncLevelDifference
  so that dPSI_mut > 0  => the inclusion isoform is MORE used in the MUTANT
  (= "inclusion_in_mut"), and dPSI_mut < 0 => "exclusion_in_mut".

THE INFERENCE BACKBONE:
  With n = 1, the high-confidence findings are genes that move in the SAME
  direction in BOTH lineages.  The EMX1 n Nestin concordant set is written out
  separately and is what steps 15-17 prioritise.

OUTPUTS  ->  results/09_AS_interpretation/genelists/
  <comp>_significant_events.csv      every event passing the filter (+ annotations)
  <comp>_genes_all.txt               unique gene symbols (any direction)
  <comp>_genes_inclusion.txt         genes with dominant inclusion_in_mut
  <comp>_genes_exclusion.txt         genes with dominant exclusion_in_mut
  <comp>_background.txt              all gene symbols rMATS could test (ORA universe)
  background_union.txt               union of both comparison backgrounds
  concordant_genes.csv               EMX1 n Nestin, same direction (high confidence)
  concordant_genes.txt               just the symbols (for enrichment / STRING)
  genelist_summary.txt               human-readable counts + sanity check

USAGE (on the cluster):
  python 13_AS_genelists.py            # uses default cluster paths below
  python 13_AS_genelists.py --base-dir /path/to/90-1265107649 --min-dpsi 0.1 --min-reads 20
"""

import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

# --- comparisons & event types in scope -------------------------------------
COMPARISONS = ["EMX1_wt_vs_mut", "Nestin_wt_vs_mut"]
EVENT_TYPES = ["SE", "A5SS", "A3SS", "RI"]   # MXE intentionally excluded

# --- curated splicing-regulator panel (for at-a-glance annotation only) ------
# Used to flag, not to filter. Drawn from the spliceosome / SR / hnRNP / neuronal
# RBP families repeatedly implicated in the project's snRNA SJU work and in
# Soelter et al. 2026 (same Setbp1 S858R model).
SPLICING_REGULATORS = {
    "Son", "Srrm1", "Srrm2", "Sf3b1", "Sf3b2", "Sf1", "Pnisr", "Srek1", "Zcchc7",
    "Hnrnpa1", "Hnrnpa2b1", "Hnrnpu", "Hnrnpk", "Hnrnpc", "Hnrnpd", "Hnrnpm",
    "Mbnl1", "Mbnl2", "Rbm25", "Rbm5", "Rbm39", "Ddx5", "Ddx17", "Gpatch8",
    "Snrpn", "Snrnp70", "Tra2a", "Tra2b", "Ptbp1", "Ptbp2", "Rbfox1", "Rbfox2",
    "Rbfox3", "Nova1", "Nova2", "Celf1", "Celf2", "Khdrbs1", "Elavl1", "Elavl2",
    "Elavl3", "Elavl4", "Srsf1", "Srsf2", "Srsf3", "Srsf5", "Srsf6", "Srsf7",
    "Srsf11", "Tut4", "Zcchc11", "Nrdc", "Nrd1",
}

# Soelter et al. 2026 cross-tissue RBPs (named in the paper text) -- a positive
# control: do our concordant hits recover these?
SOELTER_RBPS = {"Pnisr", "Hnrnpa2b1", "Srrm2", "Zcchc7", "Son", "Mbnl1", "Srek1"}


def find_count_columns(df):
    """rMATS trailing columns are identically NAMED across SE/A5SS/A3SS/RI even
    though their POSITIONS differ; selecting by name is robust to event type."""
    needed = ["IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
              "IncLevel1", "IncLevel2", "IncLevelDifference", "PValue", "FDR",
              "GeneID", "geneSymbol", "chr", "strand"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"rMATS table missing expected columns: {missing}")
    return needed


def _first_float(val):
    """IncLevel can be a comma-joined list (one per replicate); n=1 => single
    value, but parse defensively by taking the mean of present values."""
    if pd.isna(val):
        return np.nan
    parts = [p for p in str(val).split(",") if p not in ("", "NA")]
    if not parts:
        return np.nan
    try:
        return float(np.mean([float(p) for p in parts]))
    except ValueError:
        return np.nan


def load_comparison(comp_dir, min_dpsi, min_reads):
    """Return (sig_events_df, background_symbols_set) for one comparison."""
    sig_rows = []
    background = set()

    for event in EVENT_TYPES:
        path = os.path.join(comp_dir, f"{event}.MATS.JC.txt")
        if not os.path.isfile(path):
            print(f"  [warn] missing {path} -- skipping {event}", file=sys.stderr)
            continue

        df = pd.read_csv(path, sep="\t", dtype={"GeneID": str, "geneSymbol": str})
        find_count_columns(df)

        # background = every gene rMATS tested for this event type
        background.update(df["geneSymbol"].dropna().astype(str))

        for col in ["IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
                    "PValue", "FDR"]:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        inc1 = df["IncLevel1"].map(_first_float)
        inc2 = df["IncLevel2"].map(_first_float)
        ild = pd.to_numeric(df["IncLevelDifference"], errors="coerce")

        reads_s1 = df["IJC_SAMPLE_1"] + df["SJC_SAMPLE_1"]
        reads_s2 = df["IJC_SAMPLE_2"] + df["SJC_SAMPLE_2"]

        keep = (
            (ild.abs() >= min_dpsi)
            & (reads_s1 >= min_reads)
            & (reads_s2 >= min_reads)
            & inc1.notna() & inc2.notna()
        )

        sub = df.loc[keep].copy()
        if sub.empty:
            continue
        sub["event_type"] = event
        sub["IncLevel1_num"] = inc1[keep].values
        sub["IncLevel2_num"] = inc2[keep].values
        sub["dPSI_mut"] = -ild[keep].values            # >0 => inclusion in mutant
        sub["abs_dPSI"] = ild[keep].abs().values
        sub["min_reads_either"] = np.minimum(reads_s1[keep].values,
                                             reads_s2[keep].values)
        sub["direction"] = np.where(sub["dPSI_mut"] > 0,
                                    "inclusion_in_mut", "exclusion_in_mut")
        sig_rows.append(sub)

    if sig_rows:
        events = pd.concat(sig_rows, ignore_index=True)
    else:
        events = pd.DataFrame()
    return events, background


def collapse_to_genes(events):
    """One row per gene: dominant direction = sign of the largest |dPSI| event."""
    if events.empty:
        return pd.DataFrame(columns=["geneSymbol", "GeneID", "n_events",
                                     "max_abs_dPSI", "dom_dPSI_mut",
                                     "dom_direction", "is_splicing_regulator"])
    rows = []
    for sym, g in events.groupby("geneSymbol"):
        top = g.loc[g["abs_dPSI"].idxmax()]
        rows.append({
            "geneSymbol": sym,
            "GeneID": top["GeneID"],
            "n_events": len(g),
            "max_abs_dPSI": float(g["abs_dPSI"].max()),
            "dom_dPSI_mut": float(top["dPSI_mut"]),
            "dom_direction": top["direction"],
            "is_splicing_regulator": sym in SPLICING_REGULATORS,
        })
    return pd.DataFrame(rows).sort_values("max_abs_dPSI", ascending=False)


def write_lines(path, items):
    with open(path, "w") as fh:
        for it in sorted(set(items)):
            fh.write(f"{it}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-dir",
                    default="/beegfs/scratch/ric.sessa/kubacki.michal/"
                            "SRF_Linda_top/90-1265107649",
                    help="path to the 90-1265107649 project root")
    ap.add_argument("--min-dpsi", type=float, default=0.10)
    ap.add_argument("--min-reads", type=int, default=20)
    args = ap.parse_args()

    splicing_dir = os.path.join(args.base_dir, "results", "05_splicing")
    out_dir = os.path.join(args.base_dir, "results", "09_AS_interpretation",
                           "genelists")
    os.makedirs(out_dir, exist_ok=True)

    per_gene = {}          # comp -> gene-level DataFrame
    backgrounds = {}       # comp -> set
    summary = []
    summary.append("AS interpretation -- step 13 (gene lists & concordance)")
    summary.append("=" * 60)
    summary.append(f"Filter: |dPSI| >= {args.min_dpsi} AND reads >= "
                   f"{args.min_reads} in both samples (PValue/FDR NOT used).")
    summary.append("n=1 per condition: all results are EXPLORATORY; the "
                   "concordant set is the inference backbone.\n")

    for comp in COMPARISONS:
        comp_dir = os.path.join(splicing_dir, comp)
        print(f"[13] processing {comp} ...")
        events, background = load_comparison(comp_dir, args.min_dpsi,
                                             args.min_reads)
        genes = collapse_to_genes(events)
        per_gene[comp] = genes
        backgrounds[comp] = background

        # event table
        if not events.empty:
            cols = ["event_type", "GeneID", "geneSymbol", "chr", "strand",
                    "IncLevel1_num", "IncLevel2_num", "IncLevelDifference",
                    "dPSI_mut", "abs_dPSI", "direction", "min_reads_either",
                    "PValue", "FDR"]
            events[cols].sort_values("abs_dPSI", ascending=False).to_csv(
                os.path.join(out_dir, f"{comp}_significant_events.csv"),
                index=False)

        # gene lists
        write_lines(os.path.join(out_dir, f"{comp}_genes_all.txt"),
                    genes["geneSymbol"])
        write_lines(os.path.join(out_dir, f"{comp}_genes_inclusion.txt"),
                    genes.loc[genes["dom_direction"] == "inclusion_in_mut",
                              "geneSymbol"])
        write_lines(os.path.join(out_dir, f"{comp}_genes_exclusion.txt"),
                    genes.loc[genes["dom_direction"] == "exclusion_in_mut",
                              "geneSymbol"])
        write_lines(os.path.join(out_dir, f"{comp}_background.txt"), background)

        n_reg = int(genes["is_splicing_regulator"].sum())
        summary.append(f"{comp}:")
        summary.append(f"  significant events : {len(events)}")
        summary.append(f"  unique genes       : {len(genes)}")
        summary.append(f"  background (tested) : {len(background)}")
        summary.append(f"  splicing regulators among hits: {n_reg}")
        if n_reg:
            hits = genes.loc[genes['is_splicing_regulator'],
                             'geneSymbol'].tolist()
            summary.append(f"    -> {', '.join(sorted(hits))}")
        summary.append("")

    # union background
    union_bg = set().union(*backgrounds.values()) if backgrounds else set()
    write_lines(os.path.join(out_dir, "background_union.txt"), union_bg)

    # --- concordance: gene significant in BOTH, same dominant direction --------
    g0, g1 = per_gene[COMPARISONS[0]], per_gene[COMPARISONS[1]]
    if not g0.empty and not g1.empty:
        merged = g0.merge(g1, on="geneSymbol", suffixes=(
            f"_{COMPARISONS[0]}", f"_{COMPARISONS[1]}"))
        same_dir = (merged[f"dom_direction_{COMPARISONS[0]}"]
                    == merged[f"dom_direction_{COMPARISONS[1]}"])
        concordant = merged.loc[same_dir].copy()
        concordant["is_splicing_regulator"] = concordant["geneSymbol"].isin(
            SPLICING_REGULATORS)
        concordant["in_soelter_rbps"] = concordant["geneSymbol"].isin(
            SOELTER_RBPS)
        concordant["mean_max_abs_dPSI"] = concordant[
            [f"max_abs_dPSI_{COMPARISONS[0]}",
             f"max_abs_dPSI_{COMPARISONS[1]}"]].mean(axis=1)
        concordant = concordant.sort_values("mean_max_abs_dPSI",
                                             ascending=False)
        concordant.to_csv(os.path.join(out_dir, "concordant_genes.csv"),
                          index=False)
        write_lines(os.path.join(out_dir, "concordant_genes.txt"),
                    concordant["geneSymbol"])

        summary.append("CONCORDANT (EMX1 n Nestin, same direction) "
                       "= HIGH CONFIDENCE")
        summary.append("-" * 60)
        summary.append(f"  concordant genes   : {len(concordant)}")
        reg = concordant.loc[concordant["is_splicing_regulator"], "geneSymbol"]
        summary.append(f"  splicing regulators: {len(reg)}"
                       + (f" -> {', '.join(sorted(reg))}" if len(reg) else ""))
        soel = concordant.loc[concordant["in_soelter_rbps"], "geneSymbol"]
        summary.append(f"  Soelter RBPs recovered: {len(soel)}"
                       + (f" -> {', '.join(sorted(soel))}" if len(soel) else ""))
    else:
        summary.append("CONCORDANT: one comparison had no events -- skipped.")

    summary_text = "\n".join(summary)
    with open(os.path.join(out_dir, "genelist_summary.txt"), "w") as fh:
        fh.write(summary_text + "\n")
    print("\n" + summary_text)
    print(f"\n[13] outputs written to {out_dir}")


if __name__ == "__main__":
    main()
