
## Genomic Read Matching and Annotation Pipeline (GRMAP)
 
A bioinformatics pipeline for **matching sequencing reads to a reference genome and annotating the matched sequences** with genomic features.

### **📌 Overview of the Workflow**
1. **Script 1: Read Matching & Extraction**  
   - Takes **sequencing reads** (from `illumina_reads_X.fasta.gz`) and **matches them** to a reference genome (`hg38_xxx.fasta.gz`).
   - Uses **multi-threaded processing** to improve speed.
   - **Finds locations** of matches in the reference genome and **records**:
     - **Start & End positions**  
     - **Strand orientation** (Forward `F` / Reverse `R`)  
     - **Matched sequence**
     - **Upstream and downstream context (20 bp)**
   - **Outputs:** `matchedseqs.txt` → A table of matched sequences and their genomic positions.

2. **Script 2: Genomic Feature Annotation**  
   - **Uses `matchedseqs.txt`** from Script 1.
   - Loads multiple **genomic feature files**:
     - **GFF3** (Gene annotation)
     - **TSS file** (Transcription Start Sites)
     - **CpG island file** (Methylation hotspots)
     - **RepeatMasker file** (Repetitive elements like SINEs, LINEs)
   - **Annotates each matched sequence** by:
     - Finding the **closest gene** and its type.
     - Checking if the sequence is in a **CpG island**.
     - Checking if the sequence overlaps with a **repetitive element**.
   - **Outputs:** `matchedseqs_annotate.txt` → An annotated version of `matchedseqs.txt` with additional biological context.

### **🔹 Summary of the Pipeline**
| **Step** | **Script** | **Purpose** | **Output** |
|----------|-----------|-------------|------------|
| **1. Read Matching** | `match.pl` | Matches sequencing reads to the reference genome | `matchedseqs.txt` |
| **2. Genomic Annotation** | `annotate.pl` | Annotates matched sequences with genes, CpG islands, and repeat elements | `matchedseqs_annotate.txt` |

## Benchmarking Details [match.pl] (2023) 

| Read Length | Query Size | Markers Found | Runtime  | Memory Usage |
|-------------|------------|---------------|----------|--------------|
| 40          | 1,000      | 511           | 2.94 min | 334 MB       |
| 60          | 1,000      | 446           | 3.27 min | 336 MB       |
| 80          | 1,000      | 415           | 3.83 min | 337.5 MB     |
| 100         | 1,000      | 439           | 3.23 min | 341 MB       |
| 100         | 10,000     | 4,280         | 33 min   | 341 MB       | 

## Benchmarking Details [match.pl] (Jan 2025) 

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

## Author 

**Abhinav Mishra** 
mishraabhinav36@gmail.com

## License
This project is licensed under the MIT License.
