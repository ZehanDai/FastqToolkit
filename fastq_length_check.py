#!/usr/bin/env python3
"""
fastq_length_check.py — Detect abnormal-length reads in paired-end FASTQ files.

Part of FastqToolkit.

Features:
    - Automatically detect abnormal read-length thresholds (IQR or MAD).
    - Output a per-pair TSV report: read ID, R1/R2 headers, lengths, anomaly type.
    - Optionally write abnormal reads and/or normal-only reads as FASTQ.
    - Optional SQLite backend for very large files (--use-sqlite).

Design:
    1. Stream both files once to build a length distribution (very low memory).
    2. Compute low/high length thresholds per file (IQR or MAD).
    3. Collect paired records into a store (in-memory dict or SQLite).
    4. Classify each pair and write a TSV report.
    5. Optionally emit abnormal / filtered FASTQ output.
"""

import argparse
import gzip
import os
import sqlite3
import sys
from collections import Counter


# ============================================================
# I/O helpers
# ============================================================

def open_fastq(path):
    """Open plain-text or gzip-compressed FASTQ for reading."""
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')


def parse_fastq_id(header_line):
    """Return read ID (part before first space, without '@')."""
    if not header_line.startswith('@'):
        return None
    return header_line[1:].split()[0]


def base_id(full_id):
    """Drop trailing /1 or /2."""
    if full_id.endswith('/1') or full_id.endswith('/2'):
        return full_id[:-2]
    return full_id


def iter_fastq(path, error_log=None):
    """Yield (full_id, base_id, header, seq, plus, qual) for each valid record."""
    with open_fastq(path) as f:
        line_num = 0
        while True:
            header = None
            while True:
                line = f.readline()
                if not line:
                    return
                line_num += 1
                if line.startswith('@'):
                    header = line.rstrip('\n')
                    break
                if error_log:
                    error_log.write(f"SKIP_UNEXPECTED\t{path}\tline {line_num}\n")

            seq = f.readline()
            if not seq:
                if error_log:
                    error_log.write(f"TRUNCATED_NO_SEQ\t{path}\tline {line_num}\n")
                return
            seq = seq.rstrip('\n')

            plus = f.readline()
            if not plus or not plus.startswith('+'):
                if error_log:
                    error_log.write(f"TRUNCATED_NO_PLUS\t{path}\tline {line_num}\n")
                return
            plus = plus.rstrip('\n')

            qual = f.readline()
            if not qual:
                if error_log:
                    error_log.write(f"TRUNCATED_NO_QUAL\t{path}\tline {line_num}\n")
                return
            qual = qual.rstrip('\n')

            full_id = parse_fastq_id(header)
            if full_id is None:
                if error_log:
                    error_log.write(f"INVALID_HEADER\t{path}\tline {line_num}\n")
                continue

            yield full_id, base_id(full_id), header, seq, plus, qual


# ============================================================
# Step 1: length distribution
# ============================================================

def scan_length_distribution(path, error_log=None):
    """
    Stream a FASTQ file, returning (Counter: length -> count, total_reads).
    Memory: O(distinct lengths), typically < 100 entries.
    """
    dist = Counter()
    total = 0
    for _, _, _, seq, _, _ in iter_fastq(path, error_log):
        dist[len(seq)] += 1
        total += 1
    return dist, total


# ============================================================
# Step 2: threshold detection
# ============================================================

def _weighted_quantile(dist, q):
    """Compute quantile q from a Counter(value -> count)."""
    items = sorted(dist.items())
    total = sum(dist.values())
    if total == 0:
        return 0.0
    target = q * (total - 1)
    cum = 0
    for value, count in items:
        if cum + count > target:
            return float(value)
        cum += count
    return float(items[-1][0])


def compute_thresholds(dist, method='iqr', factor=1.5,
                       fixed_low=None, fixed_high=None):
    """
    Return (low, high) length thresholds.

    method='iqr'  : Q1 - factor*IQR,  Q3 + factor*IQR
    method='mad'  : median ± factor * MAD
    method='fixed': use fixed_low / fixed_high as-is
    """
    if method == 'fixed':
        if fixed_low is None or fixed_high is None:
            raise ValueError("--method fixed requires --fixed-low and --fixed-high")
        return float(fixed_low), float(fixed_high)

    if method == 'iqr':
        q1 = _weighted_quantile(dist, 0.25)
        q3 = _weighted_quantile(dist, 0.75)
        iqr = q3 - q1
        return q1 - factor * iqr, q3 + factor * iqr

    if method == 'mad':
        median = _weighted_quantile(dist, 0.5)
        dev = Counter()
        for length, count in dist.items():
            dev[abs(length - median)] += count
        mad = _weighted_quantile(dev, 0.5)
        return median - factor * mad, median + factor * mad

    raise ValueError(f"Unknown method: {method}")


# ============================================================
# Step 3: paired record collection
# ============================================================

class MemoryStore:
    """In-memory store: base_id -> {'R1': (header, length), 'R2': (header, length)}."""

    def __init__(self):
        self.data = {}

    def add(self, b_id, tag, header, length):
        self.data.setdefault(b_id, {})[tag] = (header, length)

    def iter_pairs(self):
        for b_id, rec in self.data.items():
            yield b_id, rec.get('R1'), rec.get('R2')


class SQLiteStore:
    """Disk-backed store for very large files."""

    def __init__(self, db_path):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA synchronous=OFF")
        self.conn.execute("""
            CREATE TABLE IF NOT EXISTS pairs (
                base_id   TEXT PRIMARY KEY,
                r1_header TEXT, r1_length INTEGER,
                r2_header TEXT, r2_length INTEGER
            )
        """)
        self._pending = 0

    def add(self, b_id, tag, header, length):
        self.conn.execute(
            "INSERT OR IGNORE INTO pairs (base_id) VALUES (?)", (b_id,)
        )
        col_h, col_l = ('r1_header', 'r1_length') if tag == 'R1' \
                       else ('r2_header', 'r2_length')
        self.conn.execute(
            f"UPDATE pairs SET {col_h} = ?, {col_l} = ? WHERE base_id = ?",
            (header, length, b_id),
        )
        self._pending += 1
        if self._pending >= 100_000:
            self.conn.commit()
            self._pending = 0

    def iter_pairs(self):
        cur = self.conn.execute(
            "SELECT base_id, r1_header, r1_length, r2_header, r2_length FROM pairs"
        )
        for b_id, r1h, r1l, r2h, r2l in cur:
            r1 = (r1h, r1l) if r1h is not None else None
            r2 = (r2h, r2l) if r2h is not None else None
            yield b_id, r1, r2

    def close(self):
        self.conn.commit()
        self.conn.close()


def collect_pairs(path, tag, store, error_log=None, progress_every=1_000_000):
    """Scan a FASTQ file and register each record into the store."""
    n = 0
    for _, b_id, header, seq, _, _ in iter_fastq(path, error_log):
        store.add(b_id, tag, header, len(seq))
        n += 1
        if progress_every and n % progress_every == 0:
            print(f"  [{tag}] {n:,} records indexed", file=sys.stderr)
    return n


# ============================================================
# Step 4: classification and report
# ============================================================

def classify(r1, r2, r1_low, r1_high, r2_low, r2_high):
    """
    Given (header, length) for R1 and R2 (either may be None),
    return (status, detail).
    """
    if r1 is None and r2 is None:
        return "MISSING_BOTH", "R1:missing | R2:missing"
    if r1 is None:
        return "MISSING_R1", "R1:missing | R2:present"
    if r2 is None:
        return "MISSING_R2", "R1:present | R2:missing"

    r1_bad = not (r1_low <= r1[1] <= r1_high)
    r2_bad = not (r2_low <= r2[1] <= r2_high)

    if r1_bad and r2_bad:
        return "BOTH_ABNORMAL", f"R1:{r1[1]} | R2:{r2[1]} | both out of range"
    if r1_bad:
        return "R1_ABNORMAL", f"R1:{r1[1]} | R2:{r2[1]} | R1 out of range"
    if r2_bad:
        return "R2_ABNORMAL", f"R1:{r1[1]} | R2:{r2[1]} | R2 out of range"
    return "OK", f"R1:{r1[1]} | R2:{r2[1]}"


def write_report(store, thresholds, report_path, only_abnormal=True):
    """
    Write per-pair TSV:
        common_id  R1_header  R1_length  R2_header  R2_length  status  detail
    Returns (status_counts, abnormal_ids).
    """
    r1_low, r1_high = thresholds['r1']
    r2_low, r2_high = thresholds['r2']

    status_counts = Counter()
    abnormal_ids = set()

    with open(report_path, 'w') as fo:
        fo.write("common_id\tR1_header\tR1_length\tR2_header\tR2_length\tstatus\tdetail\n")
        for b_id, r1, r2 in store.iter_pairs():
            status, detail = classify(r1, r2, r1_low, r1_high, r2_low, r2_high)
            status_counts[status] += 1
            if status != "OK":
                abnormal_ids.add(b_id)
                if only_abnormal:
                    continue
            r1h = r1[0] if r1 else "NA"
            r1l = r1[1] if r1 else "NA"
            r2h = r2[0] if r2 else "NA"
            r2l = r2[1] if r2 else "NA"
            fo.write(f"{b_id}\t{r1h}\t{r1l}\t{r2h}\t{r2l}\t{status}\t{detail}\n")

    return status_counts, abnormal_ids


# ============================================================
# Step 5: optional FASTQ output
# ============================================================

def write_subset_fastq(path, out_path, id_set, invert=False, error_log=None):
    """
    Write records whose base_id is in id_set
    (or NOT in id_set when invert=True).
    Supports .gz output if out_path ends with .gz.
    """
    is_gz = out_path.endswith('.gz')
    opener = gzip.open if is_gz else open
    mode = 'wt' if is_gz else 'w'
    n = 0
    with opener(out_path, mode) as fo:
        for _, b_id, header, seq, plus, qual in iter_fastq(path, error_log):
            hit = b_id in id_set
            if hit != invert:
                fo.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
                n += 1
    return n


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Detect abnormal-length reads in paired-end FASTQ files."
    )
    parser.add_argument("--r1", required=True, help="R1 FASTQ file (.fastq/.fastq.gz)")
    parser.add_argument("--r2", required=True, help="R2 FASTQ file (.fastq/.fastq.gz)")
    parser.add_argument("--method", choices=['iqr', 'mad', 'fixed'], default='iqr',
                        help="Threshold detection method (default: iqr)")
    parser.add_argument("--factor", type=float, default=1.5,
                        help="Multiplier for IQR or MAD (default: 1.5)")
    parser.add_argument("--fixed-low", type=int, default=None,
                        help="Lower length threshold (with --method fixed)")
    parser.add_argument("--fixed-high", type=int, default=None,
                        help="Upper length threshold (with --method fixed)")
    parser.add_argument("--report", default="length_report.tsv",
                        help="Per-pair statistics TSV (default: length_report.tsv)")
    parser.add_argument("--report-all", action="store_true",
                        help="Include OK pairs in the report (default: abnormal only)")
    parser.add_argument("--output-abnormal", default=None, metavar="PREFIX",
                        help="Write abnormal reads to PREFIX_R1.fastq / PREFIX_R2.fastq")
    parser.add_argument("--output-filtered", default=None, metavar="PREFIX",
                        help="Write normal-only reads to PREFIX_R1.fastq / PREFIX_R2.fastq")
    parser.add_argument("--use-sqlite", action="store_true",
                        help="Use SQLite backend (recommended for very large files)")
    parser.add_argument("--db", default="temp_length_check.db",
                        help="SQLite database path (only with --use-sqlite)")
    parser.add_argument("--error-log", default="length_check_errors.log",
                        help="Log file for malformed records")
    args = parser.parse_args()

    if args.method == 'fixed' and (args.fixed_low is None or args.fixed_high is None):
        print("Error: --method fixed requires --fixed-low and --fixed-high",
              file=sys.stderr)
        sys.exit(1)

    error_log = open(args.error_log, 'w')
    error_log.write("# REASON\tFILE\tDETAIL\n")

    # ---- Step 1: scan length distributions -------------------------------
    print(f"Scanning length distribution of {args.r1} ...", file=sys.stderr)
    r1_dist, r1_total = scan_length_distribution(args.r1, error_log)
    print(f"  R1 reads: {r1_total:,}, distinct lengths: {len(r1_dist)}",
          file=sys.stderr)

    print(f"Scanning length distribution of {args.r2} ...", file=sys.stderr)
    r2_dist, r2_total = scan_length_distribution(args.r2, error_log)
    print(f"  R2 reads: {r2_total:,}, distinct lengths: {len(r2_dist)}",
          file=sys.stderr)

    # ---- Step 2: compute thresholds --------------------------------------
    r1_low, r1_high = compute_thresholds(
        r1_dist, args.method, args.factor, args.fixed_low, args.fixed_high
    )
    r2_low, r2_high = compute_thresholds(
        r2_dist, args.method, args.factor, args.fixed_low, args.fixed_high
    )
    print(f"Thresholds — R1: [{r1_low:.1f}, {r1_high:.1f}]  "
          f"R2: [{r2_low:.1f}, {r2_high:.1f}]", file=sys.stderr)

    # ---- Step 3: collect paired records ----------------------------------
    if args.use_sqlite:
        if os.path.exists(args.db):
            os.remove(args.db)
        store = SQLiteStore(args.db)
    else:
        store = MemoryStore()

    print(f"Collecting R1 records ...", file=sys.stderr)
    collect_pairs(args.r1, 'R1', store, error_log)
    print(f"Collecting R2 records ...", file=sys.stderr)
    collect_pairs(args.r2, 'R2', store, error_log)

    # ---- Step 4: write report --------------------------------------------
    thresholds = {'r1': (r1_low, r1_high), 'r2': (r2_low, r2_high)}
    print(f"Writing report to {args.report} ...", file=sys.stderr)
    status_counts, abnormal_ids = write_report(
        store, thresholds, args.report, only_abnormal=not args.report_all
    )

    print("\n=== LENGTH CHECK SUMMARY ===", file=sys.stderr)
    for k, v in sorted(status_counts.items()):
        print(f"  {k}: {v:,}", file=sys.stderr)

    # ---- Step 5: optional FASTQ output -----------------------------------
    if args.output_abnormal:
        print(f"Writing abnormal reads to {args.output_abnormal}_R1.fastq ...",
              file=sys.stderr)
        write_subset_fastq(args.r1, f"{args.output_abnormal}_R1.fastq",
                           abnormal_ids, invert=False, error_log=error_log)
        write_subset_fastq(args.r2, f"{args.output_abnormal}_R2.fastq",
                           abnormal_ids, invert=False, error_log=error_log)

    if args.output_filtered:
        print(f"Writing filtered reads to {args.output_filtered}_R1.fastq ...",
              file=sys.stderr)
        write_subset_fastq(args.r1, f"{args.output_filtered}_R1.fastq",
                           abnormal_ids, invert=True, error_log=error_log)
        write_subset_fastq(args.r2, f"{args.output_filtered}_R2.fastq",
                           abnormal_ids, invert=True, error_log=error_log)

    # ---- Cleanup ---------------------------------------------------------
    if isinstance(store, SQLiteStore):
        store.close()
        if os.path.exists(args.db):
            os.remove(args.db)

    error_log.close()
    print(f"Done. Error log: {args.error_log}", file=sys.stderr)


if __name__ == "__main__":
    main()