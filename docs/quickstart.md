# Quickstart

**Clone the repo**   

```bash
git clone https://github.com/bibymaths/grmap
cd grmap
``` 
 
**Edit `config.yaml`**

- Set `read` to your reads pattern.
- Point `reference`, `annotations`, `scripts_dir`, and `results_dir` to the correct paths.

**Create your environment**

```bash
conda create -n grmap -c conda-forge -c bioconda snakemake matplotlib numpy=1.26
conda activate grmap 
sudo cpan Parallel::ForkManager IO::Uncompress::Gunzip
```
 
**Run the pipeline**

```bash
snakemake --cores 1 --report results/report.html --report-after-run
```

