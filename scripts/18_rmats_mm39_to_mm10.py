#!/usr/bin/env python3
"""
================================================================================
18_rmats_mm39_to_mm10.py  --  liftOver rMATS output from GRCm39 (mm39) to mm10
================================================================================

rMAPS2 (the RBP motif-map web server) only offers the mouse assembly **mm10**,
but this project's rMATS results are **GRCm39 = mm39**. rMAPS2 reads the
coordinates in the uploaded rMATS file and extracts sequence from *its own* mm10
genome, so mm39 coordinates would point at the wrong sequence. This script
converts the coordinate columns of each rMATS event file mm39 -> mm10 with UCSC
`liftOver`, producing upload-ready mm10 rMATS files.

WHAT IT DOES
  For each event type (SE, MXE, A5SS, A3SS, RI):
    - extracts every exon interval from the rMATS coordinate columns,
    - liftOvers them all in one pass,
    - keeps an event ONLY IF every one of its intervals lifts cleanly:
        * maps to exactly one mm10 interval (no split),
        * stays on the SAME chromosome and SAME strand,
        * preserves start < end,
    - rewrites the rMATS file with the mm10 coordinates, leaving every other
      column (counts, PSI, p-value, FDR, ...) untouched.
  Events that fail (mostly unplaced scaffolds, or build-divergent regions) are
  dropped and tallied in a manifest -- this is expected and safe; the main
  chromosomes lift at high rates.

WHY ALL-OR-NOTHING PER EVENT: an rMATS event is only meaningful if its exon and
both flanks stay in correct relative position. Lifting some coordinates but not
others would corrupt the event geometry, so we drop the whole event instead.

INPUT  : results/05_splicing/<comp>/{SE,MXE,A5SS,A3SS,RI}.MATS.JC.txt   (mm39)
OUTPUT : <out-dir>/<comp>/{SE,...}.MATS.JC.txt                          (mm10)
         <out-dir>/<comp>/liftover_manifest.csv                         (counts)
         <out-dir>/liftover_overview.txt

REQUIRES: UCSC `liftOver` on PATH (conda: bioconda::ucsc-liftover) and the
          mm39ToMm10 chain file.

USAGE (cluster):
  python 18_rmats_mm39_to_mm10.py \
      --base-dir /beegfs/.../90-1265107649 \
      --chain    /path/to/mm39ToMm10.over.chain.gz \
      --out-dir  /beegfs/.../90-1265107649/results/10_rMAPS_tracks/trackB_rMAPS2_mm10/mm10_rmats_input \
      [--liftover /path/to/liftOver]
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile

COMPARISONS = ["EMX1_wt_vs_mut", "Nestin_wt_vs_mut"]

# Coordinate column pairs (0-based start, end) per event type, as named in the
# rMATS *.MATS.JC.txt header. Each pair is one exon interval, BED-compatible.
EVENT_COORD_COLS = {
    "SE":   [("exonStart_0base", "exonEnd"), ("upstreamES", "upstreamEE"),
             ("downstreamES", "downstreamEE")],
    "MXE":  [("1stExonStart_0base", "1stExonEnd"),
             ("2ndExonStart_0base", "2ndExonEnd"),
             ("upstreamES", "upstreamEE"), ("downstreamES", "downstreamEE")],
    "A5SS": [("longExonStart_0base", "longExonEnd"), ("shortES", "shortEE"),
             ("flankingES", "flankingEE")],
    "A3SS": [("longExonStart_0base", "longExonEnd"), ("shortES", "shortEE"),
             ("flankingES", "flankingEE")],
    "RI":   [("riExonStart_0base", "riExonEnd"), ("upstreamES", "upstreamEE"),
             ("downstreamES", "downstreamEE")],
}


def run_liftover(liftover_bin, in_bed, chain, out_bed, unmapped):
    res = subprocess.run([liftover_bin, in_bed, chain, out_bed, unmapped],
                         capture_output=True, text=True)
    if res.returncode != 0 and not os.path.exists(out_bed):
        raise RuntimeError(f"liftOver failed: {res.stderr.strip()}")


def convert_event(path_in, path_out, event, liftover_bin, chain, tmpdir):
    """Returns (n_total, n_kept). Writes path_out (mm10) if any kept."""
    # Raw tab handling (NOT csv) so the original formatting -- notably the
    # double-quotes rMATS puts around GeneID/geneSymbol -- is preserved exactly.
    with open(path_in) as fh:
        lines = fh.read().splitlines()
    header_line = lines[0]
    header = header_line.split("\t")
    data = [ln.split("\t") for ln in lines[1:] if ln]
    col = {name: i for i, name in enumerate(header)}  # first occurrence wins
    # locate chr/strand and the coordinate column indices
    ci_chr, ci_strand = col["chr"], col["strand"]
    pairs = EVENT_COORD_COLS[event]
    pair_idx = [(col[s], col[e]) for s, e in pairs]

    # 1) write all intervals to a BED, name = "<rowidx>|<pairidx>"
    in_bed = os.path.join(tmpdir, f"{event}.in.bed")
    out_bed = os.path.join(tmpdir, f"{event}.out.bed")
    unmapped = os.path.join(tmpdir, f"{event}.unmapped.bed")
    with open(in_bed, "w") as bed:
        for ri, fields in enumerate(data):
            chrom, strand = fields[ci_chr], fields[ci_strand]
            for pj, (si, ei) in enumerate(pair_idx):
                try:
                    s, e = int(fields[si]), int(fields[ei])
                except (ValueError, IndexError):
                    continue
                if e <= s:
                    continue
                bed.write(f"{chrom}\t{s}\t{e}\t{ri}|{pj}\t0\t{strand}\n")

    run_liftover(liftover_bin, in_bed, chain, out_bed, unmapped)

    # 2) collect lifted intervals; detect splits (duplicate names)
    lifted, seen = {}, set()
    if os.path.exists(out_bed):
        with open(out_bed) as fh:
            for ln in fh:
                c, s, e, name, _score, strand = ln.rstrip("\n").split("\t")[:6]
                if name in seen:
                    lifted[name] = None          # split mapping -> mark bad
                    continue
                seen.add(name)
                lifted[name] = (c, int(s), int(e), strand)

    # 3) reassemble each event; keep only if ALL pairs lifted consistently
    kept = []
    n_pairs = len(pair_idx)
    for ri, fields in enumerate(data):
        orig_chr, orig_strand = fields[ci_chr], fields[ci_strand]
        ok = True
        new_coords = {}
        for pj, (si, ei) in enumerate(pair_idx):
            rec = lifted.get(f"{ri}|{pj}")
            if (not rec or rec[0] != orig_chr or rec[3] != orig_strand
                    or rec[2] <= rec[1]):
                ok = False
                break
            new_coords[(si, ei)] = (rec[1], rec[2])
        if not ok:
            continue
        row = list(fields)
        for (si, ei), (s, e) in new_coords.items():
            row[si], row[ei] = str(s), str(e)
        kept.append(row)

    if kept:
        with open(path_out, "w") as out:
            out.write(header_line + "\n")
            for row in kept:
                out.write("\t".join(row) + "\n")
    return len(data), len(kept)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-dir",
                    default="/beegfs/scratch/ric.sessa/kubacki.michal/"
                            "SRF_Linda_top/90-1265107649")
    ap.add_argument("--chain", required=True, help="mm39ToMm10.over.chain(.gz)")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--liftover", default="liftOver",
                    help="path to the UCSC liftOver binary (default: on PATH)")
    args = ap.parse_args()

    splic = os.path.join(args.base_dir, "results", "05_splicing")
    out_dir = args.out_dir or os.path.join(
        args.base_dir, "results", "10_rMAPS_tracks", "trackB_rMAPS2_mm10",
        "mm10_rmats_input")
    os.makedirs(out_dir, exist_ok=True)

    overview = ["rMAPS2 Track B -- rMATS mm39 -> mm10 liftOver",
                "=" * 60,
                f"chain: {args.chain}",
                "Per-event all-or-nothing: an event is kept only if every exon "
                "interval lifts to the same chr/strand with no split.\n"]

    with tempfile.TemporaryDirectory() as tmp:
        for comp in COMPARISONS:
            comp_out = os.path.join(out_dir, comp)
            os.makedirs(comp_out, exist_ok=True)
            manifest = [("event_type", "n_total", "n_kept", "n_dropped",
                         "pct_kept")]
            overview.append(f"{comp}:")
            for event in EVENT_COORD_COLS:
                p_in = os.path.join(splic, comp, f"{event}.MATS.JC.txt")
                if not os.path.isfile(p_in):
                    overview.append(f"  {event}: input missing -- skipped")
                    continue
                p_out = os.path.join(comp_out, f"{event}.MATS.JC.txt")
                n_tot, n_keep = convert_event(p_in, p_out, event,
                                              args.liftover, args.chain, tmp)
                pct = 100.0 * n_keep / n_tot if n_tot else 0.0
                manifest.append((event, n_tot, n_keep, n_tot - n_keep,
                                 f"{pct:.1f}"))
                overview.append(f"  {event}: {n_keep}/{n_tot} lifted "
                                f"({pct:.1f}%)")
                print(f"[18] {comp} {event}: {n_keep}/{n_tot} ({pct:.1f}%)")
            with open(os.path.join(comp_out, "liftover_manifest.csv"),
                      "w", newline="") as fh:
                csv.writer(fh).writerows(manifest)
            overview.append("")

    with open(os.path.join(out_dir, "liftover_overview.txt"), "w") as fh:
        fh.write("\n".join(overview) + "\n")
    print("\n" + "\n".join(overview))
    print(f"\n[18] mm10 rMATS files in {out_dir}")


if __name__ == "__main__":
    main()
