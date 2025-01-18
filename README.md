# Read Pattern Tool

## Overview
This Perl script matches marker sequences from Illumina reads (lengths: 40, 60, 80, 100) to the partial hg38 reference genome. It calculates how many of these markers appear in chromosome 1 of hg38. Results include matched sequences, runtime, and memory benchmarks.

## Features
- Supports `.fasta.gz` input files to save memory and time.
- Matches sequences using regex for global counts of matches.
- Handles different Illumina read lengths and query sizes.
- Outputs matched markers to `MatchedMarkers.txt`.
- Benchmarks runtime and memory usage.

## Prerequisites
### Windows
- Install Perl or ActiveState Perl.
- Optional: Install Valgrind for benchmarking.

### MacOS/Linux
- Perl is pre-installed in most distributions.
- Optional: Use `/usr/bin/time` for benchmarking.

## Installation
Clone the repository and place the following files in the same directory:
1. This script (`read_pattern.pl`)
2. Reference genome file (`hg38_partial.fasta.gz`)
3. Illumina read files:
   - `illumina_reads_40.fasta.gz`
   - `illumina_reads_60.fasta.gz`
   - `illumina_reads_80.fasta.gz`
   - `illumina_reads_100.fasta.gz`

## Usage
### Windows
Run the script using any Perl-compatible IDE or command line:
```bash
perl read_pattern.pl
```

### MacOS/Linux
#### Basic Usage:
```bash
perl read_pattern.pl
```
#### Benchmarking:
```bash
/usr/bin/time -l perl read_pattern.pl  # MacOS
/usr/bin/time -v perl read_pattern.pl  # Linux
```

## Input Options
1. Select an Illumina read file:
   - 1: `illumina_reads_40.fasta.gz`
   - 2: `illumina_reads_60.fasta.gz`
   - 3: `illumina_reads_80.fasta.gz`
   - 4: `illumina_reads_100.fasta.gz`

2. Select a query size for read length = 100:
   - 1: 1,000 queries
   - 2: 10,000 queries
   - 3: 100,000 queries
   - 4: 1,000,000 queries

## Output
- Matched marker sequences are written to `MatchedMarkers.txt`.
- The script prints the total number of matches to the console.

## Benchmarking Details
| Read Length | Query Size | Markers Found | Runtime  | Memory Usage |
|-------------|------------|---------------|----------|--------------|
| 40          | 1,000      | 511           | 2.94 min | 334 MB       |
| 60          | 1,000      | 446           | 3.27 min | 336 MB       |
| 80          | 1,000      | 415           | 3.83 min | 337.5 MB     |
| 100         | 1,000      | 439           | 3.23 min | 341 MB       |
| 100         | 10,000     | 4,280         | 33 min   | 341 MB       |

## Error Handling
- Validates `.fasta.gz` file paths.
- Prompts for re-execution if invalid options are entered.

## Author
**Abhinav Mishra** (Group 7)

## License
This project is licensed under the MIT License.

