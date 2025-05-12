
# Output Structure

After successful execution, your `results/` directory will contain, for each sample:

```
results/
├── sample1.matched       # intermediate match output (removed by default)
├── sample1.annotated     # final annotated results
├── sample2.matched
└── sample2.annotated
``` 
--- 

| File Type     | Description |
|---------------|-------------|
| **`.matched`**   | Intermediate tab-delimited file of read→reference matches. Marked `temp` in the Snakefile and deleted unless `--notemp` is used. |
| **`.annotated`** | Final annotated file per sample, with gene, TSS, CpG, and repeat annotations. Includes gene ID, feature type, TSS distance, CpG count, and repeat overlaps. |

