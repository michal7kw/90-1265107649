#!/usr/bin/env python3
"""
================================================================================
14_AS_enrichment.py  --  GO / pathway over-representation on AS gene lists
================================================================================

Over-representation analysis (ORA) of the alternatively-spliced gene lists from
step 13, mirroring the project's existing enrichment pattern
(SRF_Linda_RNA/combine_data/ANALYSIS/4_GO_Terms/raw_17_go_enrichment_final_L1.py):
gseapy -> Enrichr, organism Mouse, filter Adjusted P < 0.05, xlsx + dotplot.

TWO THINGS THIS DOES DIFFERENTLY, ON PURPOSE
  1. BACKGROUND = the rMATS-tested gene universe from step 13
     (<comp>_background.txt), NOT the whole genome.  AS detection is biased
     toward long / highly-expressed genes, so the only fair ORA universe is the
     set of genes that *could* have been called spliced.
  2. A spliceosome / RNA-processing FOCUS pass quantifies the "splicing cascade"
     signal directly: how many input genes are splicing regulators, and which
     RNA-processing terms come up.  This is the positive control for the whole
     project's hypothesis.

n = 1 CAVEAT: every term here is HYPOTHESIS-GENERATING.  Output files are tagged
accordingly.  Prioritise the `concordant` gene set (EMX1 n Nestin) -- it is the
high-confidence list.

INPUTS  (from step 13)  results/09_AS_interpretation/genelists/
  <comp>_genes_all.txt / _inclusion.txt / _exclusion.txt   + <comp>_background.txt
  concordant_genes.txt  (uses background_union.txt)

OUTPUTS ->  results/09_AS_interpretation/enrichment/
  <set>/<set>_enrichment.xlsx          all significant terms (Adj P < 0.05)
  <set>/<set>_dotplot.png              top terms
  <set>/<set>_splicing_focus.csv       RNA-processing/spliceosome terms only
  enrichment_overview.txt              per-set term counts + regulator counts

USAGE (cluster):
  python 14_AS_enrichment.py
  python 14_AS_enrichment.py --base-dir /path/to/90-1265107649
Requires internet (Enrichr API) on the submit/compute node, and `gseapy`.
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import gseapy as gp

ORGANISM = "mouse"   # gseapy >=1.x validates lowercase; "Mouse" is rejected
GENE_SETS = [
    "GO_Biological_Process_2023",
    "GO_Cellular_Component_2023",
    "GO_Molecular_Function_2023",
    "Reactome_2022",
    "KEGG_2019_Mouse",
]
PADJ = 0.05
N_TOP_PLOT = 15

# terms matching these (case-insensitive) define the "splicing cascade" focus
SPLICING_FOCUS_KEYWORDS = [
    "splic", "spliceosom", "mrna processing", "rna processing", "rna binding",
    "mrna metabolic", "snrnp", "rna splicing", "ribonucleoprotein",
    "alternative", "exon", "mrna export", "nonsense-mediated",
]

# Which step-13 gene lists to test, and which background each uses.
# (comparison label, gene-list filename, background filename)
JOBS = [
    ("EMX1_all",       "EMX1_wt_vs_mut_genes_all.txt",        "EMX1_wt_vs_mut_background.txt"),
    ("EMX1_inclusion", "EMX1_wt_vs_mut_genes_inclusion.txt",  "EMX1_wt_vs_mut_background.txt"),
    ("EMX1_exclusion", "EMX1_wt_vs_mut_genes_exclusion.txt",  "EMX1_wt_vs_mut_background.txt"),
    ("Nestin_all",       "Nestin_wt_vs_mut_genes_all.txt",       "Nestin_wt_vs_mut_background.txt"),
    ("Nestin_inclusion", "Nestin_wt_vs_mut_genes_inclusion.txt", "Nestin_wt_vs_mut_background.txt"),
    ("Nestin_exclusion", "Nestin_wt_vs_mut_genes_exclusion.txt", "Nestin_wt_vs_mut_background.txt"),
    ("concordant",     "concordant_genes.txt",                "background_union.txt"),
]


def read_list(path):
    if not os.path.isfile(path):
        return []
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def dotplot(df, title, out_png):
    """Simple dotplot of top terms, x = gene ratio, colour = -log10(adj P)."""
    top = df.sort_values("Adjusted P-value").head(N_TOP_PLOT).copy()
    if top.empty:
        return
    top["Clean_Term"] = (top["Term"].str.replace(r" \(GO:[0-9]+\)", "",
                                                  regex=True)
                         .str.replace(r" R-MMU-[0-9]+", "", regex=True)
                         .str.replace("_", " "))
    top["Gene_Count"] = top["Genes"].apply(
        lambda x: len(str(x).split(";")) if pd.notna(x) and str(x) else 0)
    top["GeneRatio"] = top["Gene_Count"] / max(top["Gene_Count"].max(), 1)
    top = top.sort_values("Adjusted P-value", ascending=False)

    plt.figure(figsize=(10, max(4, len(top) * 0.42)))
    sc = plt.scatter(top["GeneRatio"], top["Clean_Term"], s=90,
                     c=-np.log10(top["Adjusted P-value"]), cmap="Reds",
                     edgecolor="black", linewidth=0.4)
    plt.colorbar(sc, label="-log10(Adj P)")
    for _, r in top.iterrows():
        plt.text(r["GeneRatio"] + 0.01, r["Clean_Term"], int(r["Gene_Count"]),
                 va="center", fontsize=8)
    plt.xlabel("Gene ratio (count / max)")
    plt.title(title + "\n[exploratory, n=1]", fontsize=10)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-dir",
                    default="/beegfs/scratch/ric.sessa/kubacki.michal/"
                            "SRF_Linda_top/90-1265107649")
    args = ap.parse_args()

    gl_dir = os.path.join(args.base_dir, "results", "09_AS_interpretation",
                          "genelists")
    out_dir = os.path.join(args.base_dir, "results", "09_AS_interpretation",
                           "enrichment")
    os.makedirs(out_dir, exist_ok=True)

    if not os.path.isdir(gl_dir):
        sys.exit(f"[14] gene-list dir not found: {gl_dir}\n"
                 f"     Run 13_AS_genelists.py first.")

    overview = ["AS interpretation -- step 14 (GO/pathway ORA)",
                "=" * 60,
                f"Background = rMATS-tested gene universe (NOT whole genome).",
                f"Libraries: {', '.join(GENE_SETS)}",
                f"Significance: Adjusted P < {PADJ}.  EXPLORATORY (n=1).\n"]

    for label, gl_file, bg_file in JOBS:
        genes = read_list(os.path.join(gl_dir, gl_file))
        background = read_list(os.path.join(gl_dir, bg_file))
        if not genes:
            overview.append(f"{label}: no genes -- skipped.")
            continue
        # drop unmappable pure-Ensembl 'symbols' from foreground (kept in bg)
        fg = [g for g in genes if not g.startswith("ENSMUSG")]

        print(f"[14] {label}: {len(fg)} genes vs {len(background)} background")
        try:
            enr = gp.enrichr(gene_list=fg, gene_sets=GENE_SETS,
                             organism=ORGANISM, background=background or None,
                             outdir=None, cutoff=0.1, verbose=False)
        except Exception as e:                       # noqa: BLE001
            overview.append(f"{label}: enrichr FAILED ({e})")
            continue

        if not (enr and isinstance(getattr(enr, "results", None), pd.DataFrame)
                and not enr.results.empty):
            overview.append(f"{label}: no enrichment results.")
            continue

        res = enr.results.copy()
        sig = res[res["Adjusted P-value"] < PADJ].copy()
        set_dir = os.path.join(out_dir, label)
        os.makedirs(set_dir, exist_ok=True)

        sig.to_excel(os.path.join(set_dir, f"{label}_enrichment.xlsx"),
                     index=False)
        dotplot(sig, f"AS enrichment: {label}",
                os.path.join(set_dir, f"{label}_dotplot.png"))

        # --- splicing / RNA-processing focus -------------------------------
        mask = res["Term"].str.lower().apply(
            lambda t: any(k in t for k in SPLICING_FOCUS_KEYWORDS))
        focus = res[mask].sort_values("Adjusted P-value")
        focus.to_csv(os.path.join(set_dir, f"{label}_splicing_focus.csv"),
                     index=False)
        n_focus_sig = int((focus["Adjusted P-value"] < PADJ).sum())

        overview.append(
            f"{label}: {len(sig)} sig terms (Adj P<{PADJ}); "
            f"{n_focus_sig} RNA-processing/splicing terms sig "
            f"(of {len(focus)} matched).")
        if n_focus_sig:
            top_focus = focus[focus["Adjusted P-value"] < PADJ].head(3)["Term"]
            for t in top_focus:
                overview.append(f"    + {t}")

    overview_text = "\n".join(overview)
    with open(os.path.join(out_dir, "enrichment_overview.txt"), "w") as fh:
        fh.write(overview_text + "\n")
    print("\n" + overview_text)
    print(f"\n[14] outputs written to {out_dir}")


if __name__ == "__main__":
    main()
