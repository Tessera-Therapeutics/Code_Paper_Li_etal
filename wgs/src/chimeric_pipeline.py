#!/usr/bin/env python3
"""
chimeric_pipeline.py — WGS chimeric read detection and translocation junction calling

Scans a name-sorted CRAM, detects split reads, filters high-confidence
translocation junctions, and uploads results to S3.

Outputs (uploaded to <s3-base><sample>/):
  chimeric_reads.csv        — all detected chimeric events
  junction_sequences.fa     — reference context FASTA at each junction
  junction_sequences.csv    — per-junction metadata + context sequences
  junction_reads/<junc>.fa  — per-junction split read FASTA (junctions
                              with >= --min-junction-reads support only)

Example:
  python chimeric_pipeline.py \\
      --cram   GM25256-B6.name-sorted.cram \\
      --reference hg38.fasta \\
      --s3-base s3://compbio-discovery-shared/nonLTR/wgs/rw-da-004/
"""

import argparse
import logging
import os
import re
import subprocess
import sys
import tempfile
import shutil
import time
from pathlib import Path

import boto3
from urllib.parse import urlparse

import pysam
import pandas as pd


CANONICAL = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

log = logging.getLogger("chimeric")


# ── Logging setup ──────────────────────────────────────────────────────────

def setup_logging(log_file=None):
    fmt = logging.Formatter(
        fmt="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    root = logging.getLogger("chimeric")
    root.setLevel(logging.DEBUG)

    if log_file:
        # Write to file; only also add stdout handler when running in a terminal
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        root.addHandler(fh)
        if sys.stdout.isatty():
            ch = logging.StreamHandler(sys.stdout)
            ch.setLevel(logging.DEBUG)
            ch.setFormatter(fmt)
            root.addHandler(ch)
            root.info(f"Logging to {log_file}")
    else:
        # No log file: write to stdout only
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(fmt)
        root.addHandler(ch)


# ── Argument parsing ───────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Input: either a local CRAM or an S3 folder to discover CRAMs from
    g_in = p.add_mutually_exclusive_group(required=True)
    g_in.add_argument("--cram", metavar="FILE",
                      help="Local name-sorted CRAM file (single sample)")
    g_in.add_argument("--s3-cram-folder", metavar="URI",
                      help="S3 folder to discover CRAMs from "
                           "(e.g. s3://bucket/wgs/); processes all matching files")

    p.add_argument("--cram-suffix", default=".name-sorted.cram", metavar="SUFFIX",
                   help="CRAM filename suffix to match when using --s3-cram-folder "
                        "[%(default)s]")
    p.add_argument("--samples", default=None, metavar="NAME[,NAME...]",
                   help="Comma-separated list of sample names to process in batch mode; "
                        "others are skipped. Duplicates: first occurrence wins.")
    p.add_argument("--reference", required=True, metavar="FILE",
                   help="Reference FASTA (must have .fai index)")
    p.add_argument("--s3-base",   required=True, metavar="URI",
                   help="S3 base URI, e.g. s3://bucket/prefix/  "
                        "(sample subfolder is appended automatically)")
    p.add_argument("--sample", default=None, metavar="NAME",
                   help="Sample name override (single-CRAM mode only; "
                        "default: derived from filename)")
    p.add_argument("--max-reads", type=int, default=0, metavar="N",
                   help="Max CRAM records to scan; 0 = full file [%(default)s]")
    p.add_argument("--threads", type=int, default=4, metavar="N",
                   help="pysam threads [%(default)s]")
    p.add_argument("--log-file", default=None, metavar="FILE",
                   help="Write log to this file in addition to stdout "
                        "(default: <sample>_pipeline.log next to the CRAM)")
    p.add_argument("--tmp-dir", default=None, metavar="DIR",
                   help="Directory for temporary CRAM downloads in batch mode "
                        "(default: system temp; use a path with sufficient disk space)")

    g = p.add_argument_group("chimeric detection")
    g.add_argument("--min-aln-len",            type=int,   default=20,
                   help="Min read-spanning length of each alignment, bp [%(default)s]")
    g.add_argument("--max-overlap",            type=int,   default=20,
                   help="Max allowed read-coordinate overlap, bp [%(default)s]")
    g.add_argument("--max-overlap-frac",       type=float, default=0.50,
                   help="Overlap must not exceed this fraction of the shorter aln [%(default)s]")
    g.add_argument("--min-sv-size",            type=int,   default=1_000,
                   help="Ref-gap threshold: large_deletion vs small_deletion [%(default)s]")
    g.add_argument("--min-competing-aln-frac", type=float, default=0.80,
                   help="Flag read if any other aln covers >= this fraction of the read [%(default)s]")
    g.add_argument("--max-extra-alns", type=int, default=2,
                   help="Flag read if it has more than this many extra supplementary/secondary "
                        "alignments (repeat scatter signature) [%(default)s]")

    g = p.add_argument_group("unambiguous junction filter")
    g.add_argument("--min-mapq",       type=int,   default=30,
                   help="Both alignments must have MAPQ >= this [%(default)s]")
    g.add_argument("--min-flank",      type=int,   default=20,
                   help="Each alignment flank must span >= this many bp [%(default)s]")
    g.add_argument("--max-base-freq",  type=float, default=0.70,
                   help="Skip read if any single nucleotide exceeds this fraction "
                        "of the sequence (low-complexity filter) [%(default)s]")

    g = p.add_argument_group("junction calling")
    g.add_argument("--bin-size",           type=int, default=10_000,
                   help="Breakpoint binning window, bp [%(default)s]")
    g.add_argument("--min-junction-reads", type=int, default=4,
                   help="Min supporting reads to export per-junction FASTA [%(default)s]")
    g.add_argument("--context-bp",         type=int, default=30,
                   help="Reference bases to fetch on each side of breakpoint [%(default)s]")

    return p.parse_args()


# ── CRAM scanning ──────────────────────────────────────────────────────────

def to_fwd_coords(aln):
    """Return (start, end) of the aligned region in forward-read coordinates."""
    if not aln.cigartuples:
        return 0, 0
    leading_hard  = aln.cigartuples[0][1]  if aln.cigartuples[0][0]  == 5 else 0
    trailing_hard = aln.cigartuples[-1][1] if aln.cigartuples[-1][0] == 5 else 0
    rl = aln.query_length
    if aln.is_reverse:
        fwd_start = trailing_hard + (rl - aln.query_alignment_end)
        fwd_end   = trailing_hard + (rl - aln.query_alignment_start)
    else:
        fwd_start = leading_hard + aln.query_alignment_start
        fwd_end   = leading_hard + aln.query_alignment_end
    return fwd_start, fwd_end


def _is_valid_chimeric_pair(aln1, aln2, min_aln_len, max_overlap, max_overlap_frac):
    s1, e1 = to_fwd_coords(aln1)
    s2, e2 = to_fwd_coords(aln2)
    len1, len2 = e1 - s1, e2 - s2
    if len1 < min_aln_len or len2 < min_aln_len:
        return False
    overlap = min(e1, e2) - max(s1, s2)
    if overlap > max_overlap:
        return False
    if overlap > 0 and overlap > max_overlap_frac * min(len1, len2):
        return False
    return True



CANONICAL_SET = set(CANONICAL)


def _process_group(group, events, cfg):
    for is_r1 in (True, False):
        end_reads = [r for r in group if r.is_read1 == is_r1 and not r.is_secondary]
        primary   = next((r for r in end_reads
                          if not r.is_supplementary and not r.is_unmapped), None)
        if primary is None:
            continue
        # Early skip: primary must be on a canonical chrom for a translocation to matter
        if primary.reference_name not in CANONICAL_SET:
            continue
        suppls   = [r for r in end_reads if r.is_supplementary and not r.is_unmapped]
        all_alts = [r for r in group
                    if r.is_read1 == is_r1
                    and (r.is_supplementary or r.is_secondary)
                    and not r.is_unmapped]
        total_len = primary.query_length

        for suppl in suppls:
            # Skip non-canonical or same-chrom pairs immediately
            if suppl.reference_name not in CANONICAL_SET:
                continue
            if suppl.reference_name == primary.reference_name:
                continue
            if not _is_valid_chimeric_pair(primary, suppl,
                                           cfg.min_aln_len, cfg.max_overlap,
                                           cfg.max_overlap_frac):
                continue
            ps, pe = to_fwd_coords(primary)
            ss, se = to_fwd_coords(suppl)
            overlap = min(pe, se) - max(ps, ss)
            other_alts = [a for a in all_alts if a is not suppl]
            full_aln_alt = (
                len(other_alts) > cfg.max_extra_alns
                or any(
                    (to_fwd_coords(alt)[1] - to_fwd_coords(alt)[0])
                    >= cfg.min_competing_aln_frac * total_len
                    for alt in other_alts
                )
            )
            events.append({
                "read_name"        : primary.query_name,
                "read_end"         : "R1" if is_r1 else "R2",
                "chrom1"           : primary.reference_name,
                "pos1"             : primary.reference_start,
                "end1"             : primary.reference_end,
                "strand1"          : "-" if primary.is_reverse else "+",
                "mapq1"            : primary.mapping_quality,
                "read_s1"          : ps,
                "read_e1"          : pe,
                "chrom2"           : suppl.reference_name,
                "pos2"             : suppl.reference_start,
                "end2"             : suppl.reference_end,
                "strand2"          : "-" if suppl.is_reverse else "+",
                "mapq2"            : suppl.mapping_quality,
                "read_s2"          : ss,
                "read_e2"          : se,
                "overlap_bp"       : overlap,
                "event_type"       : "translocation",
                "has_full_aln_alt" : full_aln_alt,
            })


def scan_cram(cfg):
    """Scan name-sorted CRAM and return chimeric_df."""
    t0 = time.time()
    events = []
    pysam.set_verbosity(0)
    cram = pysam.AlignmentFile(cfg.cram, "rc",
                               reference_filename=cfg.reference,
                               threads=cfg.threads)
    current_name = None
    group = []
    n_reads = n_groups = 0

    for read in cram:
        if read.query_name != current_name:
            if group:
                n_groups += 1
                _process_group(group, events, cfg)
                if n_groups % 500_000 == 0:
                    log.info(f"  {n_groups/1e6:.1f}M groups | {len(events):,} chimeric "
                             f"| {time.time()-t0:.0f}s elapsed")
            current_name = read.query_name
            group = [read]
        else:
            group.append(read)
        n_reads += 1
        if cfg.max_reads and n_reads >= cfg.max_reads:
            break

    if group:
        n_groups += 1
        _process_group(group, events, cfg)
    cram.close()

    chimeric_df = pd.DataFrame(events)
    label = f"first {cfg.max_reads//1_000_000}M reads" if cfg.max_reads else "full CRAM"
    log.info(f"Scanned {n_reads:,} reads ({n_groups:,} groups) in {time.time()-t0:.0f}s  [{label}]")
    log.info(f"Chimeric pairs found: {len(chimeric_df):,}")
    return chimeric_df


# ── Pipeline steps ─────────────────────────────────────────────────────────

def filter_unambiguous(chimeric_df, cfg):
    trans_all = chimeric_df[
        (chimeric_df["event_type"] == "translocation") &
        chimeric_df["chrom1"].isin(CANONICAL) &
        chimeric_df["chrom2"].isin(CANONICAL)
    ].copy()

    ci = {c: i for i, c in enumerate(CANONICAL)}
    swap = trans_all.apply(lambda r: ci[r.chrom1] > ci[r.chrom2], axis=1)
    for ca, cb in [("chrom1","chrom2"), ("pos1","pos2"), ("end1","end2"),
                   ("strand1","strand2"), ("mapq1","mapq2"),
                   ("read_s1","read_s2"), ("read_e1","read_e2")]:
        trans_all.loc[swap, [ca, cb]] = trans_all.loc[swap, [cb, ca]].values

    trans_all["left_flank"]  = (trans_all[["read_e1","read_e2"]].min(axis=1)
                               - trans_all[["read_s1","read_s2"]].min(axis=1))
    trans_all["right_flank"] = (trans_all[["read_e1","read_e2"]].max(axis=1)
                               - trans_all[["read_s1","read_s2"]].max(axis=1))

    unambig = (
        (trans_all["mapq1"]       >= cfg.min_mapq)  &
        (trans_all["mapq2"]       >= cfg.min_mapq)  &
        (trans_all["left_flank"]  >= cfg.min_flank) &
        (trans_all["right_flank"] >= cfg.min_flank) &
        (~trans_all["has_full_aln_alt"])
    )

    log.info(f"Canonical translocation reads : {len(trans_all):,}")
    log.info(f"  Fail MAPQ>={cfg.min_mapq}         : "
             f"{(~((trans_all.mapq1>=cfg.min_mapq)&(trans_all.mapq2>=cfg.min_mapq))).sum():,}")
    log.info(f"  Fail flank>={cfg.min_flank} bp    : "
             f"{(~((trans_all.left_flank>=cfg.min_flank)&(trans_all.right_flank>=cfg.min_flank))).sum():,}")
    log.info(f"  Have competing/scattered aln : {trans_all['has_full_aln_alt'].sum():,}"
             f"  (full-read aln ≥{cfg.min_competing_aln_frac:.0%} or >{cfg.max_extra_alns} extra alns)")
    log.info(f"  Passing unambiguous filter   : {unambig.sum():,}  ({100*unambig.mean():.1f}%)")

    return trans_all[unambig].copy()


def cluster_junctions(trans, cfg):
    trans = trans.copy()
    trans["bin1"] = (trans["end1"] // cfg.bin_size) * cfg.bin_size
    trans["bin2"] = (trans["pos2"] // cfg.bin_size) * cfg.bin_size

    junctions = (
        trans
        .groupby(["chrom1", "bin1", "chrom2", "bin2"], sort=False)
        .agg(
            n_reads      = ("read_name", "nunique"),
            bp1_median   = ("end1",  "median"),
            bp1_spread   = ("end1",  lambda x: int(x.max() - x.min())),
            bp2_median   = ("pos2",  "median"),
            bp2_spread   = ("pos2",  lambda x: int(x.max() - x.min())),
            r1_support   = ("read_end", lambda x: (x == "R1").sum()),
            r2_support   = ("read_end", lambda x: (x == "R2").sum()),
            mapq1_median = ("mapq1", "median"),
            mapq2_median = ("mapq2", "median"),
        )
        .reset_index()
        .sort_values("n_reads", ascending=False)
        .reset_index(drop=True)
    )
    junctions["label"] = (
        junctions["chrom1"].str.replace("chr", "") + ":"
        + junctions["bp1_median"].astype(int).map(lambda x: f"{x:,}")
        + " → "
        + junctions["chrom2"].str.replace("chr", "") + ":"
        + junctions["bp2_median"].astype(int).map(lambda x: f"{x:,}")
    )
    junctions["tight"] = (junctions["bp1_spread"] <= 200) & (junctions["bp2_spread"] <= 200)
    log.info(f"Unique junctions: {len(junctions):,}  (tight: {junctions['tight'].sum():,})")
    return trans, junctions


# ── S3 helpers ─────────────────────────────────────────────────────────────

def _s3_cp(local, remote):
    result = subprocess.run(["aws", "s3", "cp", local, remote, "--quiet"],
                            capture_output=True, text=True)
    if result.returncode == 0:
        log.info(f"  uploaded → {remote}")
    else:
        log.error(f"  Upload error ({remote}): {result.stderr.strip()}")


def _rc(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def _is_low_complexity(seq, max_base_freq):
    if not seq:
        return True
    seq_upper = seq.upper()
    return max(seq_upper.count(b) for b in "ACGT") / len(seq_upper) > max_base_freq


def _safe_fname(label):
    return re.sub(r"[^A-Za-z0-9_\-]", "_", label.replace(" → ", "_to_").replace(",", ""))


# ── Exporters ──────────────────────────────────────────────────────────────

def export_chimeric_csv(chimeric_df, s3_prefix):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as fh:
        chimeric_df.to_csv(fh, index=False)
        tmp = fh.name
    _s3_cp(tmp, s3_prefix + "chimeric_reads.csv")
    os.unlink(tmp)
    log.info(f"  {len(chimeric_df):,} events")


def export_junction_sequences(junctions, cfg, s3_prefix):
    fasta_file = pysam.FastaFile(cfg.reference)
    rows, fasta_lines = [], []

    for _, jrow in junctions.iterrows():
        bp1, bp2 = int(jrow.bp1_median), int(jrow.bp2_median)
        try:
            seq1 = fasta_file.fetch(jrow.chrom1, max(0, bp1 - cfg.context_bp), bp1).upper()
            seq2 = fasta_file.fetch(jrow.chrom2, bp2, bp2 + cfg.context_bp).upper()
        except Exception:
            seq1 = seq2 = "N" * cfg.context_bp

        fasta_lines.append(
            f">{jrow.label.replace(' ', '_')}  n_reads={int(jrow.n_reads)}  "
            f"bp1_spread={jrow.bp1_spread}  bp2_spread={jrow.bp2_spread}  "
            f"tight={'Y' if jrow.tight else 'N'}"
        )
        fasta_lines.append(seq1 + seq2)
        rows.append({
            "label": jrow.label, "chrom1": jrow.chrom1, "bp1": bp1,
            "seq_chrom1": seq1, "chrom2": jrow.chrom2, "bp2": bp2,
            "seq_chrom2": seq2, "junction_seq": seq1 + seq2,
            "n_reads": int(jrow.n_reads), "bp1_spread": jrow.bp1_spread,
            "bp2_spread": jrow.bp2_spread, "tight": jrow.tight,
            "r1_support": int(jrow.r1_support), "r2_support": int(jrow.r2_support),
            "mapq1_median": jrow.mapq1_median, "mapq2_median": jrow.mapq2_median,
        })
    fasta_file.close()

    seq_df = pd.DataFrame(rows)
    uploads = {}
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fh:
        fh.write("\n".join(fasta_lines) + "\n")
        uploads[fh.name] = s3_prefix + "junction_sequences.fa"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as fh:
        seq_df.to_csv(fh, index=False)
        uploads[fh.name] = s3_prefix + "junction_sequences.csv"

    for local, remote in uploads.items():
        _s3_cp(local, remote)
        os.unlink(local)
    log.info(f"  {len(seq_df):,} junctions")


def export_per_junction_reads(junctions, trans, cfg, s3_prefix):
    junc_filt = junctions[junctions["n_reads"] >= cfg.min_junction_reads]
    if len(junc_filt) == 0:
        log.warning(f"No junctions with >={cfg.min_junction_reads} reads — skipping per-junction FASTAs")
        return

    log.info(f"  {len(junc_filt):,} / {len(junctions):,} junctions with "
             f">={cfg.min_junction_reads} reads")

    join_cols   = ["chrom1", "bin1", "chrom2", "bin2"]
    trans_filt  = trans.merge(junc_filt[join_cols], on=join_cols, how="inner")
    target_keys = set(zip(trans_filt["read_name"], trans_filt["read_end"] == "R1"))
    n_target = len(target_keys)
    log.info(f"  Fetching {n_target:,} read sequences from CRAM…")

    pysam.set_verbosity(0)
    sequences = {}
    n_low_complexity = 0
    t0 = time.time()
    cram = pysam.AlignmentFile(cfg.cram, "rc",
                               reference_filename=cfg.reference,
                               threads=cfg.threads)
    n_scanned = 0
    for aln in cram:
        if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
            continue
        key = (aln.query_name, aln.is_read1)
        if key not in target_keys or key in sequences:
            continue
        raw = aln.query_sequence or ""
        if _is_low_complexity(raw, cfg.max_base_freq):
            target_keys.discard(key)
            n_low_complexity += 1
            continue
        sequences[key] = _rc(raw) if aln.is_reverse else raw
        if len(sequences) == len(target_keys):
            break
        n_scanned += 1
        if n_scanned % 500_000 == 0:
            log.info(f"    scanned {n_scanned//1_000_000:.1f}M reads | "
                     f"found {len(sequences):,} / {len(target_keys):,} "
                     f"| {time.time()-t0:.0f}s")
    cram.close()
    log.info(f"  Found {len(sequences):,} / {n_target:,} sequences in {time.time()-t0:.0f}s"
             f"  ({n_low_complexity:,} dropped as low-complexity, >{cfg.max_base_freq:.0%} single base)")

    tmp_dir = Path(tempfile.mkdtemp(prefix="junc_reads_"))
    n_written = 0
    try:
        for _, jrow in junc_filt.iterrows():
            mask = (
                (trans_filt["chrom1"] == jrow.chrom1) & (trans_filt["bin1"] == jrow.bin1) &
                (trans_filt["chrom2"] == jrow.chrom2) & (trans_filt["bin2"] == jrow.bin2)
            )
            lines = []
            for _, rrow in trans_filt[mask].iterrows():
                seq = sequences.get((rrow.read_name, rrow.read_end == "R1"))
                if not seq:
                    continue
                rs1, re1 = int(rrow.read_s1), int(rrow.read_e1)
                rs2, re2 = int(rrow.read_s2), int(rrow.read_e2)
                N = len(seq)
                if rs1 > rs2:
                    seq = _rc(seq)
                    rs1, re1 = N - re1, N - rs1
                    rs2, re2 = N - re2, N - rs2
                split = (min(re1, re2) + max(rs1, rs2)) // 2
                lines.append(
                    f">{rrow.read_name}_{rrow.read_end}"
                    f"  junc={split}  {rrow.chrom1}:{rs1}-{re1}"
                    f"  {rrow.chrom2}:{rs2}-{re2}  mapq={int(rrow.mapq1)}/{int(rrow.mapq2)}"
                )
                lines.append(seq)

            if not lines:
                continue
            with open(tmp_dir / f"{_safe_fname(jrow.label)}.fa", "w") as fh:
                fh.write("\n".join(lines) + "\n")
            n_written += 1

        s3_dest = s3_prefix + "junction_reads/"
        result = subprocess.run(
            ["aws", "s3", "sync", str(tmp_dir), s3_dest, "--quiet"],
            capture_output=True, text=True,
        )
        if result.returncode == 0:
            log.info(f"  {n_written:,} FASTA files uploaded → {s3_dest}")
        else:
            log.error(f"  S3 sync error: {result.stderr.strip()}")
    finally:
        shutil.rmtree(tmp_dir)


# ── CRAM discovery ─────────────────────────────────────────────────────────

def discover_crams(s3_folder, suffix):
    """Return list of S3 URIs for CRAMs matching suffix in s3_folder."""
    parsed = urlparse(s3_folder)
    bucket = parsed.netloc
    prefix = parsed.path.lstrip("/")

    s3 = boto3.client("s3")
    paginator = s3.get_paginator("list_objects_v2")

    keys = []
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            if obj["Key"].endswith(suffix):
                keys.append(f"s3://{bucket}/{obj['Key']}")

    return sorted(keys)


def _download_cram(s3_uri, local_dir):
    """Download CRAM + .crai index to local_dir, return local CRAM path."""
    parsed = urlparse(s3_uri)
    bucket = parsed.netloc
    key    = parsed.path.lstrip("/")
    fname  = key.split("/")[-1]
    local  = Path(local_dir) / fname

    s3 = boto3.client("s3")
    log.info(f"  Downloading {fname}…")
    s3.download_file(bucket, key, str(local))

    # Try to fetch the .crai index
    for idx_key in (key + ".crai", re.sub(r"\.cram$", ".crai", key)):
        try:
            s3.head_object(Bucket=bucket, Key=idx_key)
            idx_local = str(local) + ".crai" if idx_key == key + ".crai" \
                        else str(local).replace(".cram", ".crai")
            s3.download_file(bucket, idx_key, idx_local)
            log.info(f"  Downloaded index {idx_key.split('/')[-1]}")
            break
        except Exception:
            pass

    return str(local)


# ── Main ───────────────────────────────────────────────────────────────────

def _run_one(cfg, cram_path, sample_name, s3_prefix):
    """Run the full pipeline for a single CRAM."""
    cfg.cram   = cram_path
    cfg.sample = sample_name

    log.info(f"Sample    : {cfg.sample}")
    log.info(f"CRAM      : {cfg.cram}")
    log.info(f"S3 prefix : {s3_prefix}")

    log.info("── Step 1: scan CRAM ──────────────────────────────────────────────")
    chimeric_df = scan_cram(cfg)
    if chimeric_df.empty:
        log.error("No chimeric reads found — skipping.")
        return False

    log.info("── Step 2: unambiguous filter ─────────────────────────────────────")
    trans = filter_unambiguous(chimeric_df, cfg)
    del chimeric_df
    if trans.empty:
        log.error("No reads passed the unambiguous filter — skipping.")
        return False

    log.info("── Step 3: cluster junctions ──────────────────────────────────────")
    trans, junctions = cluster_junctions(trans, cfg)

    log.info("── Step 4: upload junction sequences ──────────────────────────────")
    export_junction_sequences(junctions, cfg, s3_prefix)

    log.info("── Step 5: upload per-junction read FASTAs ────────────────────────")
    export_per_junction_reads(junctions, trans, cfg, s3_prefix)

    return True


def main():
    cfg = parse_args()
    cfg.s3_base = cfg.s3_base.rstrip("/") + "/"

    # ── Single CRAM mode ───────────────────────────────────────────────────
    if cfg.cram:
        sample = cfg.sample or re.sub(
            r"\.(name-sorted|sorted|dedup)$", "", Path(cfg.cram).stem)
        log_file = cfg.log_file or str(
            Path(cfg.cram).parent / f"{sample}_pipeline.log")
        setup_logging(log_file)

        _run_one(cfg, cfg.cram, sample, cfg.s3_base + sample + "/")
        log.info("Done.")
        return

    # ── Batch mode: discover CRAMs from S3 ────────────────────────────────
    if not cfg.log_file:
        ts = time.strftime("%Y%m%d_%H%M%S")
        cfg.log_file = f"batch_pipeline_{ts}.log"
    setup_logging(cfg.log_file)
    log.info(f"Discovering CRAMs in {cfg.s3_cram_folder} (suffix: {cfg.cram_suffix})")
    cram_uris = discover_crams(cfg.s3_cram_folder, cfg.cram_suffix)

    if not cram_uris:
        log.error(f"No CRAMs found matching *{cfg.cram_suffix} in {cfg.s3_cram_folder}")
        sys.exit(1)

    log.info(f"Found {len(cram_uris)} CRAMs:")
    for uri in cram_uris:
        log.info(f"  {uri}")

    sample_filter = (set(s.strip() for s in cfg.samples.split(","))
                     if cfg.samples else None)

    results = {}
    seen_samples = set()
    tmp_root = Path(tempfile.mkdtemp(prefix="cram_batch_", dir=cfg.tmp_dir))
    try:
        for uri in cram_uris:
            fname  = uri.split("/")[-1]
            sample = re.sub(r"\.(name-sorted|sorted|dedup)$", "",
                            re.sub(r"\.cram$", "", fname))

            if sample_filter and sample not in sample_filter:
                continue
            if sample in seen_samples:
                log.info(f"  Skipping duplicate: {sample} ({uri})")
                continue
            seen_samples.add(sample)

            i = len(seen_samples)
            n_total = len(sample_filter) if sample_filter else len(cram_uris)
            log.info(f"\n{'='*60}")
            log.info(f"[{i}/{n_total}] {sample}")
            log.info(f"{'='*60}")

            tmp_dir = tmp_root / sample
            tmp_dir.mkdir()
            try:
                local_cram = _download_cram(uri, tmp_dir)
                ok = _run_one(cfg, local_cram, sample,
                              cfg.s3_base + sample + "/")
                results[sample] = "OK" if ok else "SKIPPED"
            except Exception as e:
                log.error(f"  Failed: {e}")
                results[sample] = f"ERROR: {e}"
            finally:
                shutil.rmtree(tmp_dir)
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)

    log.info("\n── Batch summary ──────────────────────────────────────────────────")
    for sample, status in results.items():
        log.info(f"  {sample}: {status}")
    log.info("Done.")


if __name__ == "__main__":
    main()
