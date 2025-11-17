def get_negative_control_counts(wildcards):
    return(f"data/kmer_counts/{config["negative_control_sample"]}_{wildcards.gene}.txt")

rule normalize_neg:
    input:
        "data/kmer_counts/{sample}_{gene}.txt",
        get_negative_control_counts
    output:
        "data/negative_control_normalized/{sample}_{gene}.txt"
    conda:
        "../envs/tidyverse.yaml"
    params:
        negative_control=config["negative_control_sample"]
    shell:
        """
        Rscript scripts/normalize_with_negative_control.R {wildcards.gene} {wildcards.sample} {params.negative_control}
        """
