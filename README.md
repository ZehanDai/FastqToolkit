# FastqToolkit
A collection of utility scripts for processing paired-end FASTQ files.

---
## Description
FastqToolkit provides standalone Python scripts for common paired-end FASTQ tasks
- Extracting random subsets of reads while keeping R1/R2 pairing intact.
- Checking R1/R2 consistency and computing summary statistics for pair read files.
- Detecting abnormal-length reads and reporting length distributions.
- 
**Update** 2026-09-18

**Version** v0.2.0
- `fastq_subset.py`: extract random subset of paired reads, with integrated reliability statistics.
- `fastq_length_check.py`: detect abnormal-length reads, output per-pair TSV report, and generate length-distribution tables / histograms.

## Usage
### 1. Extract a random subset of paired-end reads
python fastq_subset.py \
        --r1 $fq1 \
        --r2 $fq2 \
        --percent 0.1 \
        --seed 123 \
        --out-prefix my_subset0.1 \
        --error-log mispair_ids.txt \
        --stats-output stats.txt

"--percent" sets Fraction of common reads to keep, between 0 and 1 (e.g., 0.1 = 10%)
"--seed" sets random seed for reproducible sampling
"--out-prefix" sets prefix of the output fastq file
"--error-log"	sets file name for recording malformed or mismatched records

## Testing (optional)
**Prepare test data** (require SRAtoolkit)
```
> srrID=SRR18532519
> prefetch -c $srrID
> fasterq-dump $srrID -O $srrID/
> ls ./$srrID
SRR18532519.sra  SRR18532519_1.fastq  SRR18532519_2.fastq
> cd ./$srrID
> gzip $srrID*.fastq # optional
> cd ..
```
This would download a Chlamydia SRA dataset and convert it into paired-end FASTQ files.

**Run the subsetting script**
```
> fq1=./"$srrID"/"$srrID"_1.fastq.gz
> fq2=./"$srrID"/"$srrID"_2.fastq.gz
> time python fastq_subset.py \
        --r1 $fq1 \
        --r2 $fq2 \
        --percent 0.1 \
        --seed 123 \
        --out-prefix my_subset0.1 \
        --error-log mispair_ids.txt
```

### 2. Detect abnormal-length reads and report length distributions
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
