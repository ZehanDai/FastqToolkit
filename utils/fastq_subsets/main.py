#!/usr/bin/env python3
"""
Extract a random subset of paired-end FASTQ reads.

Features:
    - Keep only common read IDs between R1 and R2.
    - Randomly sample a fraction of common IDs.
    - Tolerate malformed records and log errors.
    - Report reliability statistics for both input files.
"""

import argparse
import gzip
import random
import sys
from collections import Counter


# ---------------------------------------------------------------------------
# FASTQ parsing helpers
# ---------------------------------------------------------------------------

def parse_fastq_id(header_line):
    """Extract read identifier (before space, without '@')."""
    if not header_line.startswith('@'):
        return None
    return header_line[1:].split()[0]


def open_fastq(filepath):
    """Open plain-text or gzip-compressed FASTQ."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')


def safe_iter_fastq(filepath, error_log_handle=None, stats=None):
    """
    Generator yielding (read_id, header, seq, plus, qual) for valid records.
    Invalid records are skipped and logged.

    If `stats` is a Counter, the following keys are tallied:
        valid_records, skip_unexpected, invalid_header,
        truncated_no_seq, truncated_no_plus, invalid_plus_line,
        truncated_no_qual
    """
    def _tally(key):
        if stats is not None:
            stats[key] += 1

    with open_fastq(filepath) as f:
        line_num = 0
        while True:
            # Find next '@' line
            header = None
            while True:
                line = f.readline()
                if not line:
                    return
                line_num += 1
                if line.startswith('@'):
                    header = line.rstrip('\n')
                    break
                else:
                    if error_log_handle:
                        error_log_handle.write(
                            f"SKIP_UNEXPECTED\t{filepath}\tline {line_num}\t{line[:80].rstrip()}\n"
                        )
                    _tally('skip_unexpected')

            seq = f.readline()
            if not seq:
                if error_log_handle:
                    error_log_handle.write(
                        f"TRUNCATED_NO_SEQ\t{filepath}\tline {line_num}\tmissing sequence\n"
                    )
                _tally('truncated_no_seq')
                continue
            seq = seq.rstrip('\n')

            plus = f.readline()
            if not plus:
                if error_log_handle:
                    error_log_handle.write(
                        f"TRUNCATED_NO_PLUS\t{filepath}\tline {line_num}\tmissing '+'\n"
                    )
                _tally('truncated_no_plus')
                continue
            plus = plus.rstrip('\n')
            if not plus.startswith('+'):
                if error_log_handle:
                    error_log_handle.write(
                        f"INVALID_PLUS_LINE\t{filepath}\tline {line_num}\t{plus[:80]}\n"
                    )
                _tally('invalid_plus_line')
                continue

            qual = f.readline()
            if not qual:
                if error_log_handle:
                    error_log_handle.write(
                        f"TRUNCATED_NO_QUAL\t{filepath}\tline {line_num}\tmissing quality\n"
                    )
                _tally('truncated_no_qual')
                continue
            qual = qual.rstrip('\n')

            read_id = parse_fastq_id(header)
            if read_id is None:
                if error_log_handle:
                    error_log_handle.write(
                        f"INVALID_HEADER\t{filepath}\tline {line_num}\t{header[:80]}\n"
                    )
                _tally('invalid_header')
                continue

            _tally('valid_records')
            yield read_id, header, seq, plus, qual


def collect_ids_from_fastq(filepath, error_log, stats=None):
    """Return set of valid read IDs from a FASTQ file."""
    ids = set()
    for read_id, _, _, _, _ in safe_iter_fastq(filepath, error_log, stats=stats):
        ids.add(read_id)
    return ids


def filter_fastq(input_path, output_path, keep_ids, error_log):
    """Write records whose ID is in keep_ids to output_path."""
    with open(output_path, 'w') as fout:
        for read_id, header, seq, plus, qual in safe_iter_fastq(input_path, error_log):
            if read_id in keep_ids:
                fout.write(f"{header}\n{seq}\n{plus}\n{qual}\n")


# ---------------------------------------------------------------------------
# Statistics report
# ---------------------------------------------------------------------------

_ERROR_KEYS = [
    'skip_unexpected', 'invalid_header',
    'truncated_no_seq', 'truncated_no_plus',
    'invalid_plus_line', 'truncated_no_qual',
]


def _file_section(label, stats):
    """Per-file statistics block, as a list of lines."""
    valid = stats.get('valid_records', 0)
    invalid = sum(stats.get(k, 0) for k in _ERROR_KEYS)
    lines = [
        f"--- {label} statistics ---",
        f"  Records attempted: {valid + invalid:,}",
        f"  Valid records: {valid:,}",
        f"  Invalid records: {invalid:,}",
    ]
    for k in _ERROR_KEYS:
        if stats.get(k, 0):
            lines.append(f"    {k}: {stats[k]:,}")
    return lines


def report_stats(r1_path, r2_path, r1_stats, r2_stats, ids_r1, ids_r2,
                 stats_output=None):
    """
    Build and emit a paired-end FASTQ reliability report.

    If `stats_output` is given, the report is written to that file;
    otherwise it is printed to stdout.
    """
    common = len(ids_r1 & ids_r2)
    r1_only = len(ids_r1 - ids_r2)
    r2_only = len(ids_r2 - ids_r1)
    total = common + r1_only + r2_only

    lines = [
        "=== PAIRED-END FASTQ RELIABILITY REPORT ===",
        f"R1: {r1_path}",
        f"R2: {r2_path}",
        "",
    ]
    lines += _file_section("R1", r1_stats)
    lines.append("")
    lines += _file_section("R2", r2_stats)
    lines.append("")
    lines.append("--- Pairing statistics ---")
    lines.append(f"  Common valid IDs: {common:,}")
    lines.append(f"  R1-only valid IDs: {r1_only:,}")
    lines.append(f"  R2-only valid IDs: {r2_only:,}")
    lines.append(f"  Total valid IDs across both files: {total:,}")
    if total:
        lines.append(
            f"  Completeness (paired): {common}/{total} ({100 * common / total:.2f}%)"
        )

    report = "\n".join(lines)
    if stats_output:
        with open(stats_output, 'w') as fo:
            fo.write(report + "\n")
        print(f"Statistics written to {stats_output}", file=sys.stderr)
    else:
        print(report)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Extract random subset of paired-end FASTQ reads."
    )
    parser.add_argument("--r1", required=True, help="R1 FASTQ file (.fastq or .fastq.gz)")
    parser.add_argument("--r2", required=True, help="R2 FASTQ file (.fastq or .fastq.gz)")
    parser.add_argument("--percent", type=float, required=True,
                        help="Fraction of common reads to keep (0-1)")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed for sampling")
    parser.add_argument("--out-prefix", type=str, default="subset",
                        help="Output prefix (e.g., subset_R1.fastq)")
    parser.add_argument("--error-log", type=str, default="error_ids.log",
                        help="Log file for malformed records")
    parser.add_argument("--stats-output", type=str, default=None,
                        help="Write statistics report to this file (default: stdout)")
    args = parser.parse_args()

    if not (0 < args.percent <= 1):
        print("Error: --percent must be between 0 and 1", file=sys.stderr)
        sys.exit(1)

    if args.seed is not None:
        random.seed(args.seed)

    error_log = open(args.error_log, 'w')
    error_log.write("# REASON\tFILE\tLINE_INFO\n")

    # Collect valid IDs from both files, tracking per-file statistics.
    r1_stats = Counter()
    r2_stats = Counter()

    print(f"Scanning {args.r1} ...", file=sys.stderr)
    ids_r1 = collect_ids_from_fastq(args.r1, error_log, stats=r1_stats)
    print(f"  Valid R1 IDs: {len(ids_r1)}", file=sys.stderr)

    print(f"Scanning {args.r2} ...", file=sys.stderr)
    ids_r2 = collect_ids_from_fastq(args.r2, error_log, stats=r2_stats)
    print(f"  Valid R2 IDs: {len(ids_r2)}", file=sys.stderr)

    common = ids_r1 & ids_r2
    if not common:
        print("Error: No common valid IDs found", file=sys.stderr)
        error_log.close()
        sys.exit(1)
    print(f"Common IDs: {len(common)}", file=sys.stderr)

    # Sample
    n_keep = int(len(common) * args.percent)
    if n_keep == 0:
        print("Warning: percent too low, no reads kept", file=sys.stderr)
        error_log.close()
        sys.exit(0)
    keep_ids = set(random.sample(list(common), n_keep))
    print(f"Keeping {len(keep_ids)} reads ({args.percent * 100:.1f}%)", file=sys.stderr)

    # Write outputs
    out_r1 = f"{args.out_prefix}_R1.fastq"
    out_r2 = f"{args.out_prefix}_R2.fastq"
    print(f"Writing R1 fastq to {out_r1} ...", file=sys.stderr)
    filter_fastq(args.r1, out_r1, keep_ids, error_log)
    print(f"Writing R2 fastq to {out_r2} ...", file=sys.stderr)
    filter_fastq(args.r2, out_r2, keep_ids, error_log)

    error_log.close()

    # Statistics report (replaces checkQ_fq)
    report_stats(args.r1, args.r2, r1_stats, r2_stats, ids_r1, ids_r2,
                 args.stats_output)

    print(f"IDs of mismatch pair written to: {args.error_log}", file=sys.stderr)


if __name__ == "__main__":
    main()
