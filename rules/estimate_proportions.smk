rule estimate_proportions:
    input:
        "data/kmers_passing_threshold/{sample}_fimH.txt",
        "data/kmers_passing_threshold/{sample}_sseL.txt",
        "data/conserved_kmer_count/{sample}_fimH.txt",
        "data/conserved_kmer_count/{sample}_sseL.txt",
        "data/kmer_db/fimH_kmer_presence_absence.tsv",
        "data/kmer_db/sseL_kmer_presence_absence.tsv"
    output:
        "data/strain_proportion_estimates/{sample}.txt"
    conda:
        "../envs/tidyverse.yaml"
    shell:
        """
        Rscript scripts/estimate_serotype_frequency.R {wildcards.sample}
        """
