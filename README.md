# FastqToolkit
A collection of utility scripts for processing paired-end FASTQ files.

---
## Description
FastqToolkit provides standalone Python scripts for common paired-end FASTQ tasks
- Extracting random subsets of reads while keeping R1/R2 pairing intact.
- Checking R1/R2 consistency and computing summary statistics for pair read files.
- Detecting abnormal-length reads and reporting length distributions.
  
**Update** 2026-09-18

**Version** v0.2.0
- `fastq_subset.py`: extract random subset of paired reads, with integrated reliability statistics.
- `fastq_length_check.py`: detect abnormal-length reads, output per-pair TSV report, and generate length-distribution tables / histograms.
