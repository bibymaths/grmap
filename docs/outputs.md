
# Output Structure

After successful execution, your `results/` directory will contain, for each sample:

```
results/
├── sample1.matched       # intermediate match output (removed by default)
├── sample1.annotated     # final annotated results
├── sample2.matched
└── sample2.annotated
```

- **`.matched`**  
  - Intermediate tab-delimited file of read→reference matches  
  - Marked `temp` in the Snakefile and cleaned up on success (use `--notemp` to retain)

- **`.annotated`**  
  - Final per-sample output with appended columns from gene, TSS, CpG, and repeat annotations  
  - One line per match, with fields such as gene ID, feature type, distance to TSS, CpG count, and repeat overlap  

