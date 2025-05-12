# Configuration Parameters

Edit `config.yaml` with the following keys:

| Key           | Description                                  | Example                                         |
|---------------|----------------------------------------------|-------------------------------------------------|
| `read`        | Glob pattern for input reads                 | `data/reads/illumina_reads_*.fasta.gz`          |
| `reference`   | Reference genome FASTA                       | `data/ref/hg38_partial.fasta.gz`                |
| `annotations` | Directory containing annotation files        | `data/annotations`                              |
| `scripts_dir` | Directory containing the Perl scripts        | `scripts`                                       |
| `results_dir` | Output directory for pipeline results        | `results`                                       |
| `cores`       | Number of parallel jobs for Snakemake        | `4`                                             |

Example `config.yaml`:

```yaml
read: "data/reads/illumina_reads_*.fasta.gz"
reference: "data/ref/hg38_partial.fasta.gz"
annotations: "data/annotations"
scripts_dir: "scripts"
results_dir: "results"
cores: 4
