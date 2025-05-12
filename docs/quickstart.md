# Quickstart

1. **Clone the repo**   

```bash
git clone https://github.com/bibymaths/grmap
cd grmap
``` 

2. **Edit `config.yaml`**

- Set `read` to your reads pattern.
- Point `reference`, `annotations`, `scripts_dir`, and `results_dir` to the correct paths.

3. **Create your environment**

```bash
conda create -n grmap -c conda-forge -c bioconda snakemake matplotlib numpy=1.26
conda activate grmap
```

4. **Run the pipeline**

```bash
snakemake --cores 2
```

5. **View outputs** 

```bash
ls results/
```
