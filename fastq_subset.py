#!/usr/bin/env python3
"""
Extract a random subset of paired-end FASTQ reads.

Features:
    - Keep only common read IDs between R1 and R2.
    - Randomly sample a fraction of common IDs.
    - Tolerate malformed records and log errors.
    - Support plain-text (.fastq) and gzip-compressed (.fastq.gz) input.
"""

import argparse
import gzip
import random
import sys


def parse_fastq_id(header_line):
    """Extract read identifier (before space, without '@')."""
    if not header_line.startswith('@'):
        return None
    return header_line[1:].split()[0]


def open_fastq(filepath):
    """Open a FASTQ file, transparently handling gzip compression."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')


def safe_iter_fastq(filepath, error_log_handle=None):
    """
    Generator yielding (read_id, header_line, seq, plus, qual) for valid records.
    Invalid records are skipped and logged.
    Supports plain-text and gzip-compressed FASTQ.
    """
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
            # Read three more lines
            seq = f.readline()
            if not seq:
                if error_log_handle:
                    error_log_handle.write(
                        f"TRUNCATED\t{filepath}\tline {line_num}\tmissing sequence\n"
                    )
                continue
            seq = seq.rstrip('\n')

            plus = f.readline()
            if not plus or not plus.startswith('+'):
                if error_log_handle:
                    error_log_handle.write(
                        f"TRUNCATED\t{filepath}\tline {line_num}\tmissing/invalid '+'\n"
                    )
                continue
            plus = plus.rstrip('\n')

            qual = f.readline()
            if not qual:
                if error_log_handle:
                    error_log_handle.write(
                        f"TRUNCATED\t{filepath}\tline {line_num}\tmissing quality\n"
                    )
                continue
            qual = qual.rstrip('\n')

            read_id = parse_fastq_id(header)
            if read_id is None:
                if error_log_handle:
                    error_log_handle.write(
                        f"INVALID_HEADER\t{filepath}\tline {line_num}\t{header[:80]}\n"
                    )
                continue
            yield read_id, header, seq, plus, qual


def collect_ids_from_fastq(filepath, error_log):
    """Return set of valid read IDs from a FASTQ file."""
    ids = set()
    for read_id, _, _, _, _ in safe_iter_fastq(filepath, error_log):
        ids.add(read_id)
    return ids


def filter_fastq(input_path, output_path, keep_ids, error_log):
    """Write records with ID in keep_ids to output_path."""
    with open(output_path, 'w') as fout:
        for read_id, header, seq, plus, qual in safe_iter_fastq(input_path, error_log):
            if read_id in keep_ids:
                fout.write(f"{header}\n")
                fout.write(f"{seq}\n")
                fout.write(f"{plus}\n")
                fout.write(f"{qual}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract random subset of paired-end FASTQ reads."
    )
    parser.add_argument("--r1", required=True, help="R1 FASTQ file (.fastq or .fastq.gz)")
    parser.add_argument("--r2", required=True, help="R2 FASTQ file (.fastq or .fastq.gz)")
    parser.add_argument("--percent", type=float, required=True, help="Fraction of reads to keep (0-1)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for sampling")
    parser.add_argument("--out-prefix", type=str, default="subset", help="Output prefix (e.g., subset_R1.fastq)")
    parser.add_argument("--error-log", type=str, default="error_ids.log", help="Log file for malformed records")
    args = parser.parse_args()

    if not (0 < args.percent <= 1):
        print("Error: --percent must be between 0 and 1", file=sys.stderr)
        sys.exit(1)

    if args.seed is not None:
        random.seed(args.seed)

    error_log = open(args.error_log, 'w')
    error_log.write("# REASON\tFILE\tLINE_INFO\n")

    # Collect valid IDs from both files
    print(f"Scanning {args.r1} ...", file=sys.stderr)
    ids_r1 = collect_ids_from_fastq(args.r1, error_log)
    print(f"  Valid R1 IDs: {len(ids_r1)}", file=sys.stderr)

    print(f"Scanning {args.r2} ...", file=sys.stderr)
    ids_r2 = collect_ids_from_fastq(args.r2, error_log)
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
    print(f"Done. IDs of mismatch pair written to: {args.error_log}", file=sys.stderr)


if __name__ == "__main__":
    main()
