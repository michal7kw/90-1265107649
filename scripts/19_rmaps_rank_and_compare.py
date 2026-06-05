#!/usr/bin/env python3
"""
================================================================================
19_rmaps_rank_and_compare.py  --  rank rMAPS2 hits & find reproducible RBPs
================================================================================

Parses the rMAPS2 Motif-Map summary tables (`pVal.up.vs.bg.RNAmap.txt` and
`pVal.dn.vs.bg.RNAmap.txt`) for each lineage, ranks every RBP-motif by its
smallest positional p-value, and -- the key output -- reports the RBPs that are
significant in BOTH lineages in the SAME direction. With n=1 per condition,
cross-lineage reproducibility is the real significance filter.

rMAPS2 table layout (per lineage, per direction):
  col 1   = RBP.motif  (e.g. "HuR.TT[GT][AG]TTT"; an RBP can appear in >1 row)
  cols 2-9 = smallest motif-enrichment p-value in each of 8 metaexon windows:
    upExon-3', upExonIntron, upstreamIntron, targetExon-5', targetExon-3',
    downstreamIntron, downExonIntron, downExon-5'
  Direction (given sample1=WT, sample2=MUT, IncLevelDifference=WT-MUT):
    UP = exon more included in WT  (= skipped more in MUTANT)
    DN = exon more included in MUTANT

MULTIPLE TESTING: ~114 RBP-motif rows x 8 regions = ~912 tests per
(lineage, direction). Default genome-wide bar = Bonferroni 0.05/912 = 5.5e-5.

OUTPUT -> <results-dir>/_ranked/
  ranked_<lineage>_<dir>.csv      every RBP-motif, min p, best region, all 8 ps
  reproducible_RBPs.csv           RBP-motifs passing alpha in BOTH lineages/dir
  comparison_overview.txt         human-readable summary

USAGE:
  python 19_rmaps_rank_and_compare.py \
      --results-dir .../trackB_rMAPS2_mm10/rMAPS2-results \
      [--alpha 5.5e-5] [--lineages EMX1_wt_vs_mut NESTIN_wt_vs_mut]
"""

import argparse
import csv
import os
import sys

REGIONS = ["upExon-3'", "upExonIntron", "upstreamIntron", "targetExon-5'",
           "targetExon-3'", "downstreamIntron", "downExonIntron",
           "downExon-5'"]
DIRECTIONS = {"up": "more included in WT (skipped in MUTANT)",
              "dn": "more included in MUTANT"}


def parse_table(path):
    """Return list of dicts: {rbp_motif, rbp, ps:[8 floats], min_p, best_region}."""
    out = []
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    for line in rows[1:]:
        if not line or len(line) < 9:
            continue
        rbp_motif = line[0]
        ps = []
        for x in line[1:9]:
            try:
                ps.append(float(x))
            except ValueError:
                ps.append(1.0)
        mn = min(ps)
        out.append({
            "rbp_motif": rbp_motif,
            "rbp": rbp_motif.split(".", 1)[0],
            "ps": ps,
            "min_p": mn,
            "best_region": REGIONS[ps.index(mn)],
        })
    return sorted(out, key=lambda d: d["min_p"])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", required=True,
                    help="rMAPS2-results dir containing per-lineage subfolders")
    ap.add_argument("--lineages", nargs="*", default=None,
                    help="lineage subfolder names (default: autodetect)")
    ap.add_argument("--alpha", type=float, default=5.5e-5,
                    help="significance bar (default Bonferroni 0.05/912)")
    args = ap.parse_args()

    rdir = args.results_dir
    lineages = args.lineages or sorted(
        d for d in os.listdir(rdir)
        if os.path.isdir(os.path.join(rdir, d))
        and os.path.isfile(os.path.join(rdir, d, "pVal.up.vs.bg.RNAmap.txt")))
    if len(lineages) < 1:
        sys.exit(f"[19] no lineage subfolders with pVal tables under {rdir}")

    out_dir = os.path.join(rdir, "_ranked")
    os.makedirs(out_dir, exist_ok=True)

    # parsed[lineage][dir] = sorted list of dicts
    parsed = {lin: {} for lin in lineages}
    for lin in lineages:
        for d in ("up", "dn"):
            p = os.path.join(rdir, lin, f"pVal.{d}.vs.bg.RNAmap.txt")
            rows = parse_table(p) if os.path.isfile(p) else []
            parsed[lin][d] = rows
            # write per-lineage ranked table
            with open(os.path.join(out_dir, f"ranked_{lin}_{d}.csv"),
                      "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["rank", "rbp", "rbp_motif", "min_p", "best_region"]
                           + [f"p_{r}" for r in REGIONS])
                for i, row in enumerate(rows, 1):
                    w.writerow([i, row["rbp"], row["rbp_motif"],
                                f"{row['min_p']:.3e}", row["best_region"]]
                               + [f"{x:.3e}" for x in row["ps"]])

    # --- cross-lineage reproducibility ---------------------------------------
    # Rolled up to BASE RBP (an RBP can have several motif-variant rows): take
    # each RBP's smallest p across its motif rows, per lineage per direction.
    # This is the biologically meaningful unit. We report, per direction, every
    # RBP that beats alpha in BOTH lineages.
    def rollup(rows):
        best = {}
        for r in rows:
            cur = best.get(r["rbp"])
            if cur is None or r["min_p"] < cur["min_p"]:
                best[r["rbp"]] = r
        return best

    repro_rows = []
    if len(lineages) >= 2:
        a, b = lineages[0], lineages[1]
        for d in ("up", "dn"):
            ra, rb = rollup(parsed[a][d]), rollup(parsed[b][d])
            for rbp in set(ra) & set(rb):
                pa, pb = ra[rbp]["min_p"], rb[rbp]["min_p"]
                if pa < args.alpha and pb < args.alpha:
                    repro_rows.append({
                        "direction": d, "rbp": rbp,
                        f"min_p_{a}": pa, f"region_{a}": ra[rbp]["best_region"],
                        f"min_p_{b}": pb, f"region_{b}": rb[rbp]["best_region"],
                        "same_region": ra[rbp]["best_region"] == rb[rbp]["best_region"],
                        "max_of_the_two_p": max(pa, pb),
                    })
        repro_rows.sort(key=lambda r: (r["direction"], r["max_of_the_two_p"]))
        with open(os.path.join(out_dir, "reproducible_RBPs.csv"),
                  "w", newline="") as fh:
            if repro_rows:
                w = csv.DictWriter(fh, fieldnames=list(repro_rows[0].keys()))
                w.writeheader()
                w.writerows(repro_rows)
            else:
                fh.write("none passed alpha in both lineages\n")

    # --- overview -------------------------------------------------------------
    ov = ["rMAPS2 hit ranking & cross-lineage comparison",
          "=" * 60,
          f"lineages: {', '.join(lineages)}",
          f"alpha (Bonferroni-style): {args.alpha:.2e}",
          f"exon sets per lineage: see up/dn/bg.coord.txt counts.\n"]
    for lin in lineages:
        ov.append(f"### {lin}")
        for d in ("up", "dn"):
            rows = parsed[lin][d]
            n_sig = sum(1 for r in rows if r["min_p"] < args.alpha)
            ov.append(f"  [{d}: {DIRECTIONS[d]}]  {n_sig} RBP-motifs < alpha")
            for r in rows[:6]:
                star = "*" if r["min_p"] < args.alpha else " "
                ov.append(f"    {star} {r['rbp_motif']:<28} "
                          f"{r['min_p']:.2e}  {r['best_region']}")
        ov.append("")
    if len(lineages) >= 2:
        la, lb = lineages[0], lineages[1]
        ov.append(f"### REPRODUCIBLE -- RBP significant in BOTH lineages, "
                  f"same direction (alpha={args.alpha:.1e})")
        ov.append("(cross-lineage reproducibility = the trustworthy filter at n=1)")
        if repro_rows:
            for r in repro_rows:
                same = "SAME region" if r["same_region"] else "diff region"
                ov.append(
                    f"  [{r['direction']}] {r['rbp']:<12} "
                    f"{la.split('_')[0]} {r[f'min_p_{la}']:.1e}"
                    f"({r[f'region_{la}']}) | "
                    f"{lb.split('_')[0]} {r[f'min_p_{lb}']:.1e}"
                    f"({r[f'region_{lb}']})  [{same}]")
        else:
            ov.append("  *** NONE *** -- no RBP beats alpha in both lineages.")
            ov.append("  => the two lineages have largely DISTINCT RBP-motif "
                      "signatures (or n=1 noise prevents overlap).")

    text = "\n".join(ov)
    with open(os.path.join(out_dir, "comparison_overview.txt"), "w",
              encoding="utf-8") as fh:
        fh.write(text + "\n")
    print(text)
    print(f"\n[19] outputs -> {out_dir}")


if __name__ == "__main__":
    main()
