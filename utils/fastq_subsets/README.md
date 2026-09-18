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
