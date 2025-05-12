from snakemake.io import glob_wildcards, expand, temp

configfile: "config.yaml"

# Derive samples from reads pattern stored in config
SAMPLES = glob_wildcards(config["read"].replace("*", "{sample}")).sample

rule all:
    input:
        expand("results/{sample}.annotated", sample=SAMPLES)

rule match:
    input:
        reads=lambda wc: config["read"].replace("*", wc.sample),
        ref=config["reference"]
    output:
        temp("results/{sample}.matched")
    params:
        script=config["scripts_dir"] + "/match.pl",
        query= lambda wc: str(config["querysize"])
    shell:
        """
        perl {params.script} {input.reads} {input.ref} {params.query} > {output}
        """

rule annotate:
    input:
        matches="results/{sample}.matched",
        gff=config["annotate"] + "/hg38_chr1_geneannotation.gff3.gz",
        tss=config["annotate"] + "/hg38_chr1_tss.txt.gz",
        cpg=config["annotate"] + "/hg38_cpg.txt.gz",
        repeat=config["annotate"] + "/hg38_repeatmasker.bed.gz"
    output:
        "results/{sample}.annotated"
    params:
        script=config["scripts_dir"] + "/annotate.pl"
    shell:
        """
        perl {params.script} {input.matches} {input.gff} {input.tss} {input.cpg} {input.repeat} > {output}
        """
