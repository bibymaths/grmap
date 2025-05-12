# GRmap
  
[![Snakemake](https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

[![DOI](https://zenodo.org/badge/592387295.svg)](https://doi.org/10.5281/zenodo.15390190)

A Snakemake-driven pipeline to  

- **Match** sequencing reads against a reference  
- **Annotate** matched regions with gene, CpG, TSS, and repeat information  

**Key features**    

- Fast and efficient read matching 
- Handles compressed FASTA and GFF/BED inputs   
- Supports multiple read lengths and query sizes
- Configurable via `config.yaml`
- Outputs per-sample annotation files in `results/`

---

*Next*: [Quickstart](quickstart.md)