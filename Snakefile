from snakemake.io import glob_wildcards, expand, temp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib
matplotlib.use("Agg")

configfile: "config.yaml"
SAMPLES = glob_wildcards(config["read"].replace("*", "{sample}")).sample
RESULTS = config["results"]

rule all:
    input:
        expand(RESULTS + "/{sample}.annotated", sample=SAMPLES),
        RESULTS+ "/all_samples.annotated.tsv",
        RESULTS+ "/summary_counts.txt",
        RESULTS+ "/match_counts.png"

rule match:
    input:
        reads=lambda wc: config["read"].replace("*", wc.sample),
        ref=config["reference"]
    output:
        temp(RESULTS + "/{sample}.matched")
    params:
        script=config["scripts_dir"] + "/match.pl",
        query= lambda wc: str(config["querysize"])
    shell:
        """
        perl {params.script} {input.reads} {input.ref} {params.query} > {output}
        """

rule annotate:
    input:
        matches=RESULTS + "/{sample}.matched",
        gff=config["annotate"] + "/hg38_chr1_geneannotation.gff3.gz",
        tss=config["annotate"] + "/hg38_chr1_tss.txt.gz",
        cpg=config["annotate"] + "/hg38_cpg.txt.gz",
        repeat=config["annotate"] + "/hg38_repeatmasker.bed.gz"
    output:
        RESULTS + "/{sample}.annotated"
    params:
        script=config["scripts_dir"] + "/annotate.pl"
    shell:
        """
        perl {params.script} {input.matches} {input.gff} {input.tss} {input.cpg} {input.repeat} > {output}
        """

rule merge_annotated:
    input:
        expand(RESULTS + "/{sample}.annotated", sample=SAMPLES)
    output:
        RESULTS + "/all_samples.annotated.tsv"
    shell:
        """
        head -n1 {RESULTS}/{SAMPLES[0]}.annotated > {output}
        tail -q -n+2 {RESULTS}/*.annotated >> {output}
        """

rule summarize_matches:
    input:
        merged=RESULTS + "/all_samples.annotated.tsv"
    output:
        RESULTS + "/summary_counts.txt"
    shell:
        """
        cut -f9,10 {input.merged} \\
          | tail -n+2 \\
          | sort \\
          | uniq -c \\
          | awk '{{print $2"\t"$3"\t"$1}}' \\
          > {output}
        """

rule plot_match_counts:
    input:
        RESULTS + "/summary_counts.txt"
    output:
        RESULTS + "/match_counts.png"
    run:
        df = pd.read_csv(input[0], sep="\t", header=None,
                         names=["gene_type","gene_id","count"])
        types = df["gene_type"].unique()
        cmap = plt.get_cmap("tab10")
        color_map = {t: cmap(i % cmap.N) for i, t in enumerate(types)}
        bar_colors = df["gene_type"].map(color_map)
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.bar(df["gene_id"], df["count"], color=bar_colors)
        ax.set_xlabel("Gene ID")
        ax.set_ylabel("Count")
        ax.set_xticks(range(len(df)))
        ax.set_xticklabels(df["gene_id"], rotation=45, ha="right")
        handles = [Patch(color=color_map[t], label=t) for t in types]
        ax.legend(handles=handles, title="Gene Type")
        plt.tight_layout()
        fig.savefig(output[0], dpi=300)
        plt.close()