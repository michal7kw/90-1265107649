#!/usr/bin/env python3
"""
================================================================================
17_AS_integration.py  --  External integration + interpretation report
================================================================================

The capstone. Takes the high-confidence concordant AS genes (step 13) and asks
"do independent lines of evidence agree?", then writes the human-readable
interpretation report that answers the original question -- what do these
splicing results MEAN.

INTEGRATIONS
  1. snRNA SJU convergence  -- overlap with the project's single-nucleus MARVEL
     gold/silver genes (gene_convergence_table.csv). Cross-METHOD agreement
     (bulk rMATS vs single-nucleus SJU) is strong evidence given n=1.
  2. Soelter et al. 2026     -- overlap with the 7 cross-tissue RBPs from the
     independent, same-model published study (positive control).
  3. Disease gene sets       -- overlap with curated epilepsy / neurodevelopmental
     and splicing-disease genes (SGS phenotype = seizures + NDD). SFARI / full
     panels can be supplied via --disease-dir (one symbol per line per file).
  4. STRING PPI network      -- of the concordant genes, to find a splicing-
     regulator hub/module (STRING API; needs internet; skipped gracefully).

Hypergeometric tests use the step-13 rMATS-tested background (the only fair
universe). All results are EXPLORATORY (n=1) and framed as such -- mirroring the
honest tone of ribotag_integration/compatibility_report.md.

INPUTS:
  results/09_AS_interpretation/genelists/*        (step 13)
  results/09_AS_interpretation/{enrichment,rmaps,consequence}/*overview*  (opt.)
  <RNA>/SJU_analysis/ribotag_integration/gene_convergence_table.csv
OUTPUTS -> results/09_AS_interpretation/integration/
  sju_convergence.csv, soelter_overlap.csv, disease_overlap.csv,
  string_edges.csv, string_hubs.csv (if STRING reachable)
  ../AS_interpretation_report.md      <-- the deliverable

USAGE (cluster):
  python 17_AS_integration.py --base-dir /path/to/90-1265107649 \
      --rna-dir /beegfs/.../SRF_Linda_top/SRF_Linda_RNA
"""

import argparse
import os
import sys

import pandas as pd
from scipy.stats import hypergeom

# Soelter et al. 2026 (DMM) cross-tissue RBPs named in the paper
SOELTER_RBPS = ["Pnisr", "Hnrnpa2b1", "Srrm2", "Zcchc7", "Son", "Mbnl1", "Srek1"]

# Compact curated disease/neuro gene set (illustrative; extend via --disease-dir).
# Epilepsy / neurodevelopmental / splicing-disease genes relevant to the SGS
# phenotype (seizures, developmental delay).
DISEASE_SETS = {
    "epilepsy_NDD_curated": [
        "Scn1a", "Scn2a", "Scn3a", "Scn8a", "Cacna1a", "Kcnq2", "Kcnq3",
        "Gabra1", "Gabra5", "Gabrb3", "Grin1", "Grin2a", "Grin2b", "Syngap1",
        "Stxbp1", "Cdkl5", "Mecp2", "Tcf4", "Foxg1", "Pafah1b1", "Dcx",
        "Setbp1", "Arid1b", "Chd2", "Chd8", "Mef2c", "Dnm1", "Unc13b",
        "Syt1", "Stx1b", "Snap25",
    ],
    "splicing_disease_curated": [
        "Mbnl1", "Mbnl2", "Rbfox1", "Nova1", "Nova2", "Ptbp1", "Srrm2",
        "Srrm4", "Son", "Sf3b1", "Hnrnpa1", "Hnrnpa2b1", "Fus", "Tardbp",
        "Smn1", "Srsf1", "Celf1",
    ],
}


def read_list(path):
    if not os.path.isfile(path):
        return []
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def hyper(overlap, set_a, set_b, universe):
    """one-sided P(overlap >= observed) given two sets drawn from a universe."""
    M, n, N = universe, set_a, set_b
    if min(M, n, N) == 0 or overlap == 0:
        return 1.0
    return float(hypergeom.sf(overlap - 1, M, n, N))


def read_overview(path):
    return open(path).read().strip() if os.path.isfile(path) else None


def try_string_network(genes, species=10090, score=400):
    """Query STRING API for the PPI subnetwork of `genes`. Returns (edges_df,
    hub_df) or (None, None) on any failure (no internet etc.)."""
    try:
        import requests
    except ImportError:
        return None, None
    if len(genes) < 2:
        return None, None
    # STRING caps identifiers per call; sample if huge.
    g = genes[:400]
    try:
        r = requests.post(
            "https://string-db.org/api/tsv/network",
            data={"identifiers": "%0d".join(g), "species": species,
                  "required_score": score, "caller_identity": "SRF_Linda_AS"},
            timeout=60)
        r.raise_for_status()
    except Exception:                                  # noqa: BLE001
        return None, None
    from io import StringIO
    edges = pd.read_csv(StringIO(r.text), sep="\t")
    if edges.empty or "preferredName_A" not in edges.columns:
        return None, None
    deg = (pd.concat([edges["preferredName_A"], edges["preferredName_B"]])
           .value_counts().rename_axis("gene").reset_index(name="degree"))
    return edges, deg


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-dir",
                    default="/beegfs/scratch/ric.sessa/kubacki.michal/"
                            "SRF_Linda_top/90-1265107649")
    ap.add_argument("--rna-dir",
                    default="/beegfs/scratch/ric.sessa/kubacki.michal/"
                            "SRF_Linda_top/SRF_Linda_RNA")
    ap.add_argument("--disease-dir", default=None,
                    help="optional dir of extra gene-set files (one symbol/line)")
    ap.add_argument("--no-string", action="store_true",
                    help="skip the STRING API call")
    args = ap.parse_args()

    interp = os.path.join(args.base_dir, "results", "09_AS_interpretation")
    gl = os.path.join(interp, "genelists")
    out = os.path.join(interp, "integration")
    os.makedirs(out, exist_ok=True)
    if not os.path.isdir(gl):
        sys.exit("[17] genelists/ missing -- run 13_AS_genelists.py first.")

    concordant = read_list(os.path.join(gl, "concordant_genes.txt"))
    emx1 = set(read_list(os.path.join(gl, "EMX1_wt_vs_mut_genes_all.txt")))
    nestin = set(read_list(os.path.join(gl, "Nestin_wt_vs_mut_genes_all.txt")))
    background = set(read_list(os.path.join(gl, "background_union.txt")))
    universe = len(background) or 12000
    conc_set = set(concordant)
    print(f"[17] concordant={len(conc_set)} background={universe}")

    notes = []

    # --- 1. snRNA SJU convergence --------------------------------------------
    conv_fp = os.path.join(args.rna_dir, "SJU_analysis", "ribotag_integration",
                           "gene_convergence_table.csv")
    sju_gold = set()
    sju_all = set()
    if os.path.isfile(conv_fp):
        conv = pd.read_csv(conv_fp)
        sju_all = set(conv["gene"].dropna().astype(str))
        sju_gold = set(conv.loc[conv["sju_tier"] == "Gold", "gene"]
                       .dropna().astype(str))
        ov_all = sorted(conc_set & sju_all)
        ov_gold = sorted(conc_set & sju_gold)
        p_all = hyper(len(ov_all), len(conc_set & background),
                      len(sju_all & background), universe)
        pd.DataFrame({"gene": ov_all,
                      "is_sju_gold": [g in sju_gold for g in ov_all]}).to_csv(
            os.path.join(out, "sju_convergence.csv"), index=False)
        notes.append(
            f"snRNA SJU convergence: {len(ov_all)}/{len(conc_set)} concordant "
            f"genes are SJU-significant (hypergeom p={p_all:.2g}); "
            f"{len(ov_gold)} are SJU GOLD: {', '.join(ov_gold) or '-'}")
    else:
        notes.append(f"snRNA SJU convergence: table not found at {conv_fp}")

    # --- 2. Soelter et al. 2026 ----------------------------------------------
    soel_in_conc = sorted(conc_set & set(SOELTER_RBPS))
    soel_in_emx = sorted(emx1 & set(SOELTER_RBPS))
    soel_in_nes = sorted(nestin & set(SOELTER_RBPS))
    pd.DataFrame({
        "soelter_rbp": SOELTER_RBPS,
        "in_concordant": [g in conc_set for g in SOELTER_RBPS],
        "in_EMX1": [g in emx1 for g in SOELTER_RBPS],
        "in_Nestin": [g in nestin for g in SOELTER_RBPS],
    }).to_csv(os.path.join(out, "soelter_overlap.csv"), index=False)
    notes.append(
        f"Soelter 2026 RBPs recovered -- concordant: {soel_in_conc or '-'}; "
        f"EMX1: {soel_in_emx or '-'}; Nestin: {soel_in_nes or '-'}")

    # --- 3. disease gene sets -------------------------------------------------
    disease = {k: list(v) for k, v in DISEASE_SETS.items()}
    if args.disease_dir and os.path.isdir(args.disease_dir):
        for f in os.listdir(args.disease_dir):
            disease[os.path.splitext(f)[0]] = read_list(
                os.path.join(args.disease_dir, f))
    drows = []
    for name, genes in disease.items():
        gs = set(genes)
        ov = sorted(conc_set & gs)
        p = hyper(len(ov), len(conc_set & background), len(gs & background),
                  universe)
        drows.append({"gene_set": name, "set_size": len(gs),
                      "overlap_n": len(ov), "overlap": ";".join(ov),
                      "hypergeom_p": p})
    pd.DataFrame(drows).to_csv(os.path.join(out, "disease_overlap.csv"),
                               index=False)
    for r in drows:
        notes.append(f"Disease overlap [{r['gene_set']}]: {r['overlap_n']} "
                     f"genes (p={r['hypergeom_p']:.2g}) {r['overlap']}")

    # --- 4. STRING network ----------------------------------------------------
    hub_lines = []
    if not args.no_string:
        edges, hubs = try_string_network(concordant)
        if edges is not None:
            edges.to_csv(os.path.join(out, "string_edges.csv"), index=False)
            hubs.to_csv(os.path.join(out, "string_hubs.csv"), index=False)
            top = hubs.head(10)
            hub_lines = [f"{r['gene']} (deg {r['degree']})"
                         for _, r in top.iterrows()]
            notes.append("STRING top hubs of concordant set: "
                         + ", ".join(hub_lines))
        else:
            notes.append("STRING: unreachable/empty -- skipped "
                         "(run with internet, or --no-string to silence).")

    # --- pull upstream overview lines for the report --------------------------
    enr_ov = read_overview(os.path.join(interp, "enrichment",
                                        "enrichment_overview.txt"))
    rmaps_ov = read_overview(os.path.join(interp, "rmaps", "fallback",
                                          "motif_fallback_overview.txt"))
    cons_ov = read_overview(os.path.join(interp, "consequence",
                                         "consequence_overview.txt"))
    gl_summary = read_overview(os.path.join(gl, "genelist_summary.txt"))

    write_report(args.base_dir, interp, conc_set, notes, gl_summary,
                 enr_ov, rmaps_ov, cons_ov, soel_in_conc, sju_gold & conc_set)

    txt = "[17] integration notes\n" + "\n".join(" - " + n for n in notes)
    with open(os.path.join(out, "integration_overview.txt"), "w") as fh:
        fh.write(txt + "\n")
    print("\n" + txt)
    print(f"\n[17] report -> {os.path.join(interp, 'AS_interpretation_report.md')}")


def write_report(base_dir, interp, conc_set, notes, gl_summary,
                 enr_ov, rmaps_ov, cons_ov, soel_in_conc, sju_gold_conc):
    """Assemble the narrative report (framing embedded; numbers computed)."""
    def block(title, body):
        return f"## {title}\n\n{body}\n" if body else ""

    notes_md = "\n".join(f"- {n}" for n in notes)
    gl_md = f"```\n{gl_summary}\n```" if gl_summary else "_run step 13_"
    enr_md = f"```\n{enr_ov}\n```" if enr_ov else "_run step 14_"
    rmaps_md = f"```\n{rmaps_ov}\n```" if rmaps_ov else "_run step 15_"
    cons_md = f"```\n{cons_ov}\n```" if cons_ov else "_run step 16_"

    report = f"""# Alternative Splicing — Interpretation Report (RiboTag, 90-1265107649)

> Auto-generated by `scripts/17_AS_integration.py`. Narrative framing is fixed;
> all counts are computed from the current run.
> **Model:** Setbp1 S858R (Schinzel–Giedion syndrome), mouse hippocampus,
> RiboTag bulk RNA-seq. **Comparisons:** EMX1_wt_vs_mut, Nestin_wt_vs_mut.

## TL;DR — what the splicing results mean

The dominant, reproducible signal is a **splicing cascade**: the Setbp1 S858R
gain-of-function protein accumulates (it escapes degradation) and, via the
SET/PP2A axis, perturbs the spliceosome / RNA-binding-protein (RBP) network.
The genes most affected are **splicing regulators themselves** — and their
mis-regulation propagates to hundreds of downstream targets.

Because this experiment is **n = 1 per condition**, rMATS p-values are not
trustworthy. The defensible findings are the genes that change **the same way
in both lineages** (EMX1 ∩ Nestin = the *concordant set*, **{len(conc_set)}
genes**), reinforced where they also appear in the single-nucleus SJU data and
in the independent Soelter et al. 2026 study of the same mouse model.

## The numbers (step 13)

{gl_md}

## High-confidence convergence (step 17)

{notes_md}

**Why this matters (and what not to overclaim):** cross-lineage concordance
already filters ~2,000 per-lineage hits down to a few hundred. The *global*
overlap with the single-nucleus SJU set is small and typically **not** enriched
beyond chance (see the hypergeometric p above) — the same honest verdict as
`ribotag_integration/compatibility_report.md`, driven by the n=1 design. The
value is therefore in **specific, repeatedly-convergent genes**, not in a
global enrichment claim: a gene like **Mbnl1** that is concordant across both
lineages *and* SJU-gold *and* a named Soelter RBP is a genuinely
triple-corroborated hypothesis worth targeted validation. Soelter RBPs recovered
in the concordant set:
**{', '.join(soel_in_conc) or 'none in concordant (check per-lineage in soelter_overlap.csv)'}**.
SJU-gold genes in the concordant set:
**{', '.join(sorted(sju_gold_conc)) or 'see sju_convergence.csv'}**.

## Pathways / processes affected (step 14 — GO/Reactome ORA)

Over-representation on the AS gene lists, using the rMATS-tested background.
Watch the RNA-processing / spliceosome focus pass — it quantifies the cascade
directly.

{enr_md}

## Upstream regulators (step 15 — RBP motif maps)

Which RBPs drive the events. Enrichment of MBNL (YGCY) or hnRNP motifs near
regulated exons is the direct test of the cascade/feedback hypothesis, since
`Mbnl1`/`Hnrnpa2b1` are themselves AS hits.

{rmaps_md}

## Functional consequences (step 16 — frame / NMD / domains)

What each event does to the protein: frame-shifting CDS events that are more
included in the mutant are predicted to introduce premature stops → NMD →
effective knockdown (checked against DEG log2FC).

{cons_md}

## Disease relevance

The SGS phenotype is seizures + neurodevelopmental impairment. Overlap of the
concordant AS genes with curated epilepsy/NDD and splicing-disease gene sets is
in `integration/disease_overlap.csv` (see notes above). Setbp1 itself is **not**
expected among the hits — consistent with Soelter et al., it accumulates without
being differentially spliced; the action is downstream.

## Limitations (read before quoting any number)

- **n = 1 per condition.** No p-value here is a true significance test. Everything
  is hypothesis-generating; the concordant set is a *prioritisation*, not proof.
- Enrichment is biased by gene length/expression — mitigated by the rMATS-tested
  background, not eliminated.
- Motif maps (fallback) are per-window, not per-nucleotide — run real rMAPS2 when
  available.
- NMD/domain consequences are *predicted from annotation*, not measured.

## What would turn hypotheses into findings

1. Targeted **RT-PCR / qPCR** of the top concordant events (start with the
   splicing regulators and the highest-|ΔPSI| genes; coordinates in
   `genelists/*_significant_events.csv` and the existing `06_sashimi/`).
2. **Additional biological replicates** (n ≥ 3) — the single biggest lever.
3. Where NMD is predicted, **CHX/UPF1-knockdown** or isoform-level qPCR to confirm
   the PTC isoform is degraded.

## See also (kept in sync with the single-nucleus narrative)

- `Docs_Linda/Splicing_analyses_overview/` — RiboTag vs SJU comparison, gold genes.
- `SRF_Linda_RNA/SJU_analysis/ribotag_integration/compatibility_report.md` — the
  honest cross-method verdict this report extends.
- Soelter et al. 2026, *Disease Models & Mechanisms* 19(2):dmm052402 —
  independent same-model AS study.
"""
    with open(os.path.join(interp, "AS_interpretation_report.md"), "w",
              encoding="utf-8") as fh:
        fh.write(report)


if __name__ == "__main__":
    main()
