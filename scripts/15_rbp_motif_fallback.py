#!/usr/bin/env python3
"""
================================================================================
15_rbp_motif_fallback.py  --  RBP motif enrichment around regulated exons
================================================================================

A self-contained stand-in for rMAPS2 for when the rMAPS2 standalone is not
installed on the cluster.  It answers the same biological question rMAPS2 does,
more crudely: are known RNA-binding-protein binding motifs ENRICHED in the
intronic regions flanking the exons that change splicing in the mutant, versus
the regions flanking all exons rMATS tested?

Scope: SE (skipped-exon) events for the two mutation comparisons -- SE is the
dominant, best-resolved event type and the one with the cleanest flanking-intron
geometry.

METHOD
  Foreground = SE events passing the step-13 filter (|dPSI|>=0.1, reads>=20).
  Background = all SE events rMATS tested for that comparison.
  For each event we extract four strand-corrected windows:
      - upstream intron, WIN nt immediately 5' of the cassette exon's 3' SS
      - downstream intron, WIN nt immediately 3' of the cassette exon's 5' SS
      - the cassette exon body (capped at WIN)
  For each RBP motif (degenerate IUPAC -> regex) we compute the fraction of
  windows containing >=1 occurrence in foreground vs background, and a one-sided
  Fisher exact test (foreground-enriched).  This is a *map-lite*: positional
  resolution is coarse (per-window, not per-nucleotide), so treat it as a
  hypothesis screen, subordinate to a real rMAPS2 run when available.

  Because Mbnl1 / Hnrnpa2b1 are themselves AS hits, enrichment of MBNL (YGCY)
  or hnRNP motifs here is the direct test of the splicing-cascade / feedback
  hypothesis.

DEPENDENCIES: pyfaidx (pip install pyfaidx), scipy, pandas, numpy, matplotlib.
GENOME: GRCm39 primary assembly FASTA (same assembly as the rMATS GTF).

OUTPUTS -> results/09_AS_interpretation/rmaps/fallback/
  <comp>_motif_enrichment.csv     per motif per region: fg/bg fractions, OR, p
  <comp>_motif_barplot.png        -log10(p) per motif/region
  motif_fallback_overview.txt

USAGE (cluster):
  python 15_rbp_motif_fallback.py --base-dir /path/to/90-1265107649 \
      --genome /beegfs/.../refdata-gex-GRCm39-2024-A/fasta/genome.fa
"""

import argparse
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

COMPARISONS = ["EMX1_wt_vs_mut", "Nestin_wt_vs_mut"]
WIN = 250           # nt window in each flanking intron / capped exon length
MIN_DPSI = 0.10
MIN_READS = 20

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
         "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
         "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}

# RBP motifs as DNA (transcript-sense). Families repeatedly implicated in the
# project's snRNA SJU work + Soelter et al. (Mbnl1, hnRNPs) are first.
RBP_MOTIFS = {
    "MBNL_YGCY":      "YGCY",       # Mbnl1/2 -- itself an AS hit here
    "RBFOX_UGCAUG":   "TGCATG",     # Rbfox1/2/3
    "NOVA_YCAY":      "YCAY",       # Nova1/2
    "PTBP_UCUU":      "TCTT",       # Ptbp1/2 (CU-rich)
    "PTBP_CUrich":    "YCYCY",      # Ptbp polypyrimidine
    "HNRNPA1_UAGG":   "TAGG",       # hnRNP A1
    "HNRNPA1_GGGG":   "GGGG",       # G-run, hnRNP A1/F/H
    "CELF_UGUU":      "TGTT",       # Celf1/2 (CUGBP)
    "QKI_ACUAAY":     "ACTAAY",     # Quaking
    "SRSF1_GGAGGA":   "GGAGGA",     # SR ESE
    "TIA_Urich":      "TTTTT",      # Tia1/Tial1 U-rich
}


def comp_re(motif):
    return re.compile("".join(IUPAC[b] for b in motif))


COMPILED = {name: comp_re(m) for name, m in RBP_MOTIFS.items()}
RC = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def revcomp(s):
    return s.translate(RC)[::-1]


def _first_float(val):
    if pd.isna(val):
        return np.nan
    parts = [p for p in str(val).split(",") if p not in ("", "NA")]
    try:
        return float(np.mean([float(p) for p in parts])) if parts else np.nan
    except ValueError:
        return np.nan


def load_se(comp_dir):
    path = os.path.join(comp_dir, "SE.MATS.JC.txt")
    df = pd.read_csv(path, sep="\t")
    for c in ["IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
              "exonStart_0base", "exonEnd", "upstreamEE", "downstreamES"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    ild = pd.to_numeric(df["IncLevelDifference"], errors="coerce")
    r1 = df["IJC_SAMPLE_1"] + df["SJC_SAMPLE_1"]
    r2 = df["IJC_SAMPLE_2"] + df["SJC_SAMPLE_2"]
    df["is_fg"] = ((ild.abs() >= MIN_DPSI) & (r1 >= MIN_READS)
                   & (r2 >= MIN_READS))
    return df.dropna(subset=["exonStart_0base", "exonEnd",
                             "upstreamEE", "downstreamES"])


def windows(row, genome):
    """Return strand-corrected (upstream_intron, downstream_intron, exon) seqs."""
    chrom = str(row["chr"])
    if chrom not in genome:
        return None
    strand = row["strand"]
    es, ee = int(row["exonStart_0base"]), int(row["exonEnd"])
    up_ee, dn_es = int(row["upstreamEE"]), int(row["downstreamES"])
    clen = len(genome[chrom])

    def grab(a, b):
        a, b = max(0, a), min(clen, b)
        return str(genome[chrom][a:b]).upper() if b > a else ""

    # genomic-coordinate windows (5'->3' on + strand)
    up_intron = grab(max(up_ee, es - WIN), es)          # intron before exon
    dn_intron = grab(ee, min(dn_es, ee + WIN))          # intron after exon
    exon = grab(es, min(ee, es + WIN))
    if strand == "-":
        # on minus strand transcript-sense flips up/down and reverse-comps
        up_intron, dn_intron = revcomp(dn_intron), revcomp(up_intron)
        exon = revcomp(exon)
    return up_intron, dn_intron, exon


def count_hits(seqs_by_region, motif_re):
    """seqs_by_region: dict region->list[str]; return dict region->#hit windows."""
    out = {}
    for region, seqs in seqs_by_region.items():
        out[region] = sum(1 for s in seqs if s and motif_re.search(s))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-dir",
                    default="/beegfs/scratch/ric.sessa/kubacki.michal/"
                            "SRF_Linda_top/90-1265107649")
    ap.add_argument("--genome", required=True,
                    help="GRCm39 primary assembly FASTA (matching the rMATS GTF)")
    args = ap.parse_args()

    try:
        from pyfaidx import Fasta
    except ImportError:
        sys.exit("[15-fallback] pyfaidx not installed: pip install pyfaidx")

    genome = Fasta(args.genome, sequence_always_upper=True)
    splic = os.path.join(args.base_dir, "results", "05_splicing")
    out_dir = os.path.join(args.base_dir, "results", "09_AS_interpretation",
                           "rmaps", "fallback")
    os.makedirs(out_dir, exist_ok=True)

    overview = ["AS interpretation -- step 15 fallback (RBP motif enrichment)",
                "=" * 60,
                "map-lite: per-window motif presence, fg(regulated) vs bg(tested).",
                "Subordinate to a real rMAPS2 run; hypothesis screen only.\n"]
    regions = ["upstream_intron", "downstream_intron", "exon"]

    for comp in COMPARISONS:
        print(f"[15-fallback] {comp} ...")
        df = load_se(os.path.join(splic, comp))
        fg = df[df["is_fg"]]
        bg = df  # background = all tested SE events (incl. fg)

        def seqs_for(frame):
            acc = {r: [] for r in regions}
            for _, row in frame.iterrows():
                w = windows(row, genome)
                if w is None:
                    continue
                for r, s in zip(regions, w):
                    acc[r].append(s)
            return acc

        fg_seqs, bg_seqs = seqs_for(fg), seqs_for(bg)
        n_fg = {r: len(v) for r, v in fg_seqs.items()}
        n_bg = {r: len(v) for r, v in bg_seqs.items()}

        rows = []
        for name, rgx in COMPILED.items():
            fg_hit = count_hits(fg_seqs, rgx)
            bg_hit = count_hits(bg_seqs, rgx)
            for r in regions:
                a, b = fg_hit[r], n_fg[r] - fg_hit[r]
                c, d = bg_hit[r], n_bg[r] - bg_hit[r]
                if min(n_fg[r], n_bg[r]) == 0:
                    continue
                try:
                    orr, p = fisher_exact([[a, b], [c, d]],
                                          alternative="greater")
                except ValueError:
                    orr, p = np.nan, 1.0
                rows.append({
                    "motif": name, "iupac": RBP_MOTIFS[name], "region": r,
                    "fg_frac": a / max(n_fg[r], 1), "bg_frac": c / max(n_bg[r], 1),
                    "fg_n": n_fg[r], "bg_n": n_bg[r],
                    "odds_ratio": orr, "p_fisher_greater": p,
                })
        res = pd.DataFrame(rows)
        res["enrichment"] = res["fg_frac"] / res["bg_frac"].replace(0, np.nan)
        res = res.sort_values("p_fisher_greater")
        res.to_csv(os.path.join(out_dir, f"{comp}_motif_enrichment.csv"),
                   index=False)

        # barplot of -log10 p, grouped by motif/region
        plot = res.copy()
        plot["nlp"] = -np.log10(plot["p_fisher_greater"].clip(lower=1e-300))
        plot["label"] = plot["motif"] + " | " + plot["region"]
        plot = plot.sort_values("nlp").tail(20)
        plt.figure(figsize=(9, max(4, len(plot) * 0.35)))
        plt.barh(plot["label"], plot["nlp"], color="steelblue")
        plt.axvline(-np.log10(0.05), color="red", ls="--", lw=1,
                    label="p=0.05")
        plt.xlabel("-log10(Fisher p, fg-enriched)")
        plt.title(f"RBP motif enrichment (fallback): {comp}\n[hypothesis screen]",
                  fontsize=10)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"{comp}_motif_barplot.png"), dpi=150)
        plt.close()

        sig = res[res["p_fisher_greater"] < 0.05]
        overview.append(f"{comp}: {len(fg)} fg / {len(bg)} bg SE events; "
                        f"{len(sig)} motif-region pairs enriched (p<0.05).")
        for _, r in sig.head(5).iterrows():
            overview.append(f"    {r['motif']} @ {r['region']}: "
                            f"fg {r['fg_frac']:.2f} vs bg {r['bg_frac']:.2f}, "
                            f"p={r['p_fisher_greater']:.1e}")

    txt = "\n".join(overview)
    with open(os.path.join(out_dir, "motif_fallback_overview.txt"), "w") as fh:
        fh.write(txt + "\n")
    print("\n" + txt)
    print(f"\n[15-fallback] outputs in {out_dir}")


if __name__ == "__main__":
    main()
