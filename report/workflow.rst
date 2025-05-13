GRmap: Genome Read‐Mapping & Annotation Pipeline
================================================

This report documents the execution of GRmap, a Snakemake workflow for
matching short reads against a reference genome and richly annotating the
matches with gene, TSS, CpG island, and repeat information.

Configuration
-------------

Below are the settings used for this run (from `config.yaml`):

.. csv-table::
   :header: "Parameter", "Value"
   :widths: 25, 75

   "Read pattern",        "{{ snakemake.config['read'] }}"
   "Reference genome",    "{{ snakemake.config['reference'] }}"
   "Query size",          "{{ snakemake.config['querysize'] }}"
   "Annotation dir",      "{{ snakemake.config['annotate'] }}"
   "Scripts dir",         "{{ snakemake.config['scripts_dir'] }}"
   "Results dir",         "{{ snakemake.config['results'] }}"

Execution Steps
---------------

The workflow consists of the following main rules:

1. **match**
   Uses `match.pl` to scan each sample’s reads against the reference.
   Parallelized over available CPU cores.

2. **annotate**
   Runs `annotate.pl` to add gene, TSS, CpG, and repeat annotations
   to the matched sequences.

3. **merge_annotated**
   Concatenates all per‐sample annotation files into
   `all_samples.annotated.tsv`.

4. **summarize_matches**
   Produces `summary_counts.txt` with counts per gene and CpG GC content.

5. **plot_counts_and_cpg_gc**
   Generates two bar plots: gene match counts and nonzero CpG GC contents.

6. **plot_tss_distance**
   Histograms distances of matches to their nearest TSS, by gene type.

7. **plot_tss_type**
   Boxplots of TSS distances stratified by gene type.

8. **multiqc**
   Aggregates logs and QC metrics across all samples into
   the `multiqc_report.html`.

How to Reproduce
----------------

Run the entire workflow and build this report with:

.. code-block:: bash

    snakemake --cores 4 --report report/report.html

This will execute all rules, generate the figures and tables under
`{{ snakemake.config['results'] }}`, and assemble them into a single
HTML report at `report/report.html`.
