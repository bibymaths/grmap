All Samples Annotation Summary
==============================

This file (`{{ snakemake.output[0] }}`) aggregates all per‐sample `.annotated` files
into a single TSV.  It preserves the header from the first file and then
concatenates all subsequent files, stripping their headers.
