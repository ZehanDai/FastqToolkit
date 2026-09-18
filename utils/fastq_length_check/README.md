### Detect abnormal-length reads and report length distributions
Detects abnormal read lengths using an automatic threshold (IQR or MAD), reports per-pair statistics, and optionally writes length-distribution tables and histograms for both R1 and R2.

```
python fastq_length_check.py \
    --r1 $fq1 \
    --r2 $fq2 \
    --report length_report.tsv \
    --length-report len_dist \
    --length-top 20 \
    --length-png \
    --output-abnormal abnormal \
    --output-filtered clean \
    --error-log length_check_errors.log
```
| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `--r1` | yes | – | R1 FASTQ file (`.fastq` or `.fastq.gz`) |
| `--r2` | yes | – | R2 FASTQ file (`.fastq` or `.fastq.gz`) |
| `--method` | no | `iqr` | Threshold method: `iqr`, `mad`, or `fixed` |
| `--factor` | no | `1.5` | Multiplier for IQR or MAD |
| `--fixed-low` | no | `None` | Lower length threshold (with `--method fixed`) |
| `--fixed-high` | no | `None` | Upper length threshold (with `--method fixed`) |
| `--report` | no | `length_report.tsv` | Per-pair statistics TSV |
| `--report-all` | no | off | Include `OK` pairs in the report (default: abnormal only) |
| `--length-report` | no | `None` | Write length distribution to `<PREFIX>.tsv` and `<PREFIX>.txt` |
| `--length-top` | no | `None` | Show only the top N most frequent lengths in the histogram (TSV keeps all) |
| `--length-png` | no | off | Also write `<PREFIX>.png` (requires `matplotlib`) |
| `--output-abnormal` | no | `None` | Write abnormal reads to `<PREFIX>_R1.fastq` / `<PREFIX>_R2.fastq` |
| `--output-filtered` | no | `None` | Write normal-only reads to `<PREFIX>_R1.fastq` / `<PREFIX>_R2.fastq` |
| `--use-sqlite` | no | off | Use SQLite backend (recommended for very large files) |
| `--db` | no | `temp_length_check.db` | SQLite database path (with `--use-sqlite`) |
| `--error-log` | no | `length_check_errors.log` | Log file for malformed records |


Example output of length distribution
```text
length  R1_count  R1_percent  R2_count  R2_percent
149     12        0.0012      15        0.0015
150     999988    99.9988     999985    99.9985
```
Or
```
=== LENGTH DISTRIBUTION ===

--- R1 length distribution ---
  length         count    percent  R1
     149            12    0.0012%  █
     150       999,988   99.9988%  ████████████████████████████████████████

--- R2 length distribution ---
  length         count    percent  R2
     149            15    0.0015%  █
     150       999,985   99.9985%  ████████████████████████████████████████
```
