
from matplotlib.ticker import MaxNLocator
from snakemake.io import glob_wildcards, expand, temp, report, directory
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib
matplotlib.use("Agg")

report: "report/workflow.rst"
configfile: "config.yaml"
SAMPLES = glob_wildcards(config["read"].replace("*", "{sample}")).sample
RESULTS = config["results"]

rule all:
    input:
        expand(RESULTS + "/{sample}.annotated", sample=SAMPLES),
        RESULTS+ "/all_samples.annotated.tsv",
        RESULTS+ "/summary_counts.txt",
        RESULTS+ "/gene_counts.png",
        RESULTS+ "/gene_cpg_gc.png",
        RESULTS+ "/tss_distance.png",
        RESULTS+ "/tss_type.png"

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
        report(
            RESULTS + "/{sample}.annotated",
            caption="report/annotate_sample.rst",
            category="Annotation",
            subcategory="{sample}"
        )
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
        report(
            RESULTS + "/all_samples.annotated.tsv",
            caption="report/annotated_summary.rst",
            category="Summary"
        )
    shell:
        """
        head -n1 {RESULTS}/{SAMPLES[0]}.annotated > {output}
        tail -q -n+2 {RESULTS}/*.annotated >> {output}
        """

rule summarize_matches:
    input:
        merged=RESULTS + "/all_samples.annotated.tsv"
    output:
        report(
            RESULTS + "/summary_counts.txt",
            caption="report/summary_counts.rst",
            category="Summary"
        )
    shell:
        """
        cut -f9,10,15 {input.merged} \\
          | tail -n+2 \\
          | sort \\
          | uniq -c \\
          | awk '{{print $2"\t"$3"\t"$4"\t"$1}}' \\
          > {output}
        """

rule plot_counts_and_cpg_gc:
    input:
        RESULTS + "/summary_counts.txt"
    output:
        report(
            RESULTS + "/gene_counts.png",
            caption="report/gene_counts.rst",
            category="Plots",
            labels={"measure": "count"}
        ),
        report(
            RESULTS + "/gene_cpg_gc.png",
            caption="report/gene_cpg_gc.rst",
            category="Plots",
            labels={"measure": "cpg_gc_content"}
        )
    run:
        df = pd.read_csv(
            input[0],sep="\t",header=None,
            names=["gene_type", "gene_id", "cpg_gc_content", "count"]
        )

        types = df["gene_type"].unique()
        cmap = plt.get_cmap("tab10")
        color_map = {t: cmap(i % cmap.N) for i, t in enumerate(types)}
        bar_colors = df["gene_type"].map(color_map)

        fig, ax = plt.subplots(figsize=(16, 8))

        df = df.sort_values(by="count", ascending=False)

        ax.bar(df["gene_id"], df["count"], color=bar_colors)

        ax.set_xlabel(
            "Gene",
            fontsize=10,
            fontstyle="normal",
            fontfamily="serif"
        )

        ax.set_ylabel(
            "Count",
            fontsize=10,
            fontstyle="normal",
            fontfamily="serif"
        )

        if df.shape[0] < 50:
            ax.tick_params(
                axis="x",
                labelrotation=45,
                labelsize=8,
                labelcolor="black"
            )
            plt.setp(ax.get_xticklabels(),ha="right")
        else:
            ax.tick_params(axis="x",labelbottom=False)

        ax.tick_params(
            axis="y",
            labelsize=6,
            labelcolor="black"
        )
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        handles = [Patch(color=color_map[t], label=t) for t in types]
        legend = ax.legend(
            handles=handles,
            title="Type",
            fontsize=6,
            title_fontsize=10,
            frameon=True,
            loc="best"
        )

        for text in legend.get_texts():
            text.set_fontfamily("sans-serif")
            text.set_fontstyle("normal")

        plt.tight_layout()
        fig.savefig(output[0], dpi=300)
        plt.close()

        nz = df[df["cpg_gc_content"] != 0]
        if not nz.empty:
            fig2, ax2 = plt.subplots(figsize=(8, 8))
            bar_colors2 = nz["gene_type"].map(color_map)
            ax2.bar(nz["gene_id"], nz["cpg_gc_content"], color=bar_colors2)

            ax2.set_xlabel("Gene", fontsize=10, fontstyle="normal", fontfamily="serif")
            ax2.set_ylabel("CpG GC Content (%)", fontsize=10, fontstyle="normal", fontfamily="serif")
            ax2.tick_params(axis="x", labelrotation=45, labelsize=8, labelcolor="black")
            ax2.tick_params(axis="y", labelsize=6, labelcolor="black")

            legend2 = ax2.legend(handles=handles,
                                 title="Type",
                                 fontsize=6,
                                 title_fontsize=10,
                                 frameon=True,
                                 loc="best")
            for txt in legend2.get_texts():
                txt.set_fontfamily("sans-serif")

            plt.tight_layout()
            fig2.savefig(output[1], dpi=300)
            plt.close(fig2)
        else:
            # Ensure file exists to not break downstream rules
            with open(output[1],"w") as f:
                f.write("# No data available for CpG GC content\n")

rule plot_tss_distance:
    input:
        RESULTS + "/all_samples.annotated.tsv"
    output:
        report(
            RESULTS + "/tss_distance.png",
            caption="report/tss_distance.rst",
            category="Plots",
            labels={"plot": "histogram"}
        )
    run:
        df = pd.read_csv(input[0], sep="\t", comment="#", header=0)
        df = df[(df["TSS_Distance"] != "N/A") & (df["TSS_Gene_Type"] != "N/A")]
        df["TSS_Distance"] = pd.to_numeric(df["TSS_Distance"])
        df = df[df["TSS_Distance"] <= 10000]

        plt.figure(figsize=(10, 10))
        gene_types = df["TSS_Gene_Type"].unique()
        cmap = plt.get_cmap("tab10")
        color_map = {gt: cmap(i % cmap.N) for i, gt in enumerate(gene_types)}

        for gt in gene_types:
            subset = df[df["TSS_Gene_Type"] == gt]
            plt.hist(subset["TSS_Distance"], bins=50, alpha=0.6,
                     label=gt, color=color_map[gt], edgecolor="black")

        plt.xlabel("Distance (bp)", fontsize=10)
        plt.ylabel("Number of Matches", fontsize=10)
        plt.title("Proximity of Matches to TSS", fontsize=12)
        plt.grid(axis='y', linestyle=':', alpha=0.2)
        plt.legend(title="Type", fontsize=8, title_fontsize=10)
        plt.tight_layout()
        plt.savefig(output[0], dpi=300)
        plt.close()

rule plot_tss_type:
    input:
        RESULTS + "/all_samples.annotated.tsv"
    output:
        report(
            RESULTS + "/tss_type.png",
            caption="report/tss_type.rst",
            category="Plots",
            labels={"plot": "boxplot"}
        )
    run:
        df = pd.read_csv(input[0], sep="\t", comment="#", header=0)
        df = df[(df["TSS_Distance"] != "N/A") & (df["TSS_Gene_Type"] != "N/A")]
        df["TSS_Distance"] = pd.to_numeric(df["TSS_Distance"])
        df = df[df["TSS_Distance"] <= 10000]
        plt.figure(figsize=(10, 10))
        sns.boxplot(
            data=df,
            x="TSS_Gene_Type",
            y="TSS_Distance",
            hue="TSS_Gene_Type",
            palette="tab10",
            legend=False
        )
        plt.xticks(rotation=45)
        plt.ylabel("Distance (bp)")
        plt.xlabel("Type")
        plt.title("Transcription Start Site Distances", fontsize=12)
        plt.tight_layout()
        plt.savefig(output[0], dpi=300)
        plt.close()

