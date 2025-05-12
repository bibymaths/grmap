# Process Details

## 1. Matching Reads

**Script:** `scripts/match.pl`   

**Inputs:**

| Input       | Description                                               |
|-------------|-----------------------------------------------------------|
| `reads`     | Sample FASTA (e.g. `data/reads/illumina_reads_<length>.fasta.gz`) |
| `reference` | Genome FASTA (e.g. `data/ref/hg38_partial.fasta.gz`)      |


**Output:**  
- `results/{sample}.matched` (tab-delimited list of matched regions)

**Description:**  
It takes each input read set and searches against the reference, emitting one line per successful match with coordinates and basic alignment metrics.

--- 

### Benchmarking  

#### [match.pl] (Nov 2023) 

| Read Length | Query Size | Markers Found | Runtime  | Memory Usage |
|-------------|------------|---------------|----------|--------------|
| 40          | 1,000      | 511           | 2.94 min | 334 MB       |
| 60          | 1,000      | 446           | 3.27 min | 336 MB       |
| 80          | 1,000      | 415           | 3.83 min | 337.5 MB     |
| 100         | 1,000      | 439           | 3.23 min | 341 MB       |
| 100         | 10,000     | 4,280         | 33 min   | 341 MB       | 

#### [match.pl] (Jan 2025) 

| Read Length | Query Size | Markers Found | Runtime  | Memory Usage |
|-------------|------------|---------------|----------|--------------|
| 40 | 1000 | 2590 | 0.03 min | 198.16 MB |
| 40 | 10000 | 39694 | 0.23 min | 202.75 MB |
| 40 | 100000 | 400556 | 2.18 min | 248.42 MB |
| 40 | 1000000 | 400556 | 2.22 min | 248.42 MB |
| 60 | 1000 | 1395 | 0.03 min | 248.42 MB |
| 60 | 10000 | 12669 | 0.23 min | 248.42 MB |
| 60 | 100000 | 129502 | 2.12 min | 250.07 MB |
| 60 | 1000000 | 129502 | 2.10 min | 250.07 MB |
| 80 | 1000 | 1092 | 0.03 min | 250.07 MB |
| 80 | 10000 | 9866 | 0.22 min | 250.07 MB |
| 80 | 100000 | 101426 | 2.11 min | 253.08 MB |
| 80 | 1000000 | 101426 | 2.10 min | 253.08 MB |
| 100 | 1000 | 867 | 0.03 min | 253.08 MB |
| 100 | 10000 | 8822 | 0.23 min | 253.08 MB |
| 100 | 100000 | 90293 | 2.14 min | 256.41 MB |
| 100 | 1000000 | 90293 | 2.13 min | 256.41 MB |

--- 

## 2. Annotating Matches

**Script:** `scripts/annotate.pl`   

**Inputs:**

| File | Description |
|------|-------------|
| `results/{sample}.matched` | Matched read output from previous step |
| `hg38_chr1_geneannotation.gff3.gz` | Gene features |
| `hg38_chr1_tss.txt.gz` | Transcription start sites |
| `hg38_cpg.txt.gz` | CpG locations |
| `hg38_repeatmasker.bed.gz` | Repetitive elements |


**Output:**  
- `results/{sample}.annotated` (tab-delimited, with annotation columns)

**Description:**  
It enriches each match record by querying overlap with gene models, TSS positions, CpG sites, and repeat regions.  
The output adds columns for: 

- Gene ID and feature type  
- Distance to nearest TSS  
- CpG count within the matched interval  
- Overlap with repeat elements  
