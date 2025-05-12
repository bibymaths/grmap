# GRmap

A Snakemake-driven pipeline to
1. **Match** sequencing reads against a reference  
2. **Annotate** matched regions with gene, CpG, TSS, and repeat information  

**Key features**  
- Modular: separate match and annotate steps  
- Configurable via `config.yaml`  
- Handles compressed FASTA and GFF/BED inputs  
- Outputs per-sample annotation files in `results/`

---

*Next*: [Quickstart](quickstart.md)