rule kounta_input:
    output:
        "data/kmer_db/{gene}_input_alleles.txt"
    params:
        allele_path=config["allele_path"]
    shell:
        """
        for f in {params.allele_path}/{wildcards.gene}/*.fasta; do echo $f; done > {output}
        """

rule kounta:
    input:
        "data/kmer_db/{gene}_input_alleles.txt"
    output:
        kmers="data/kmer_db/{gene}_kounta.tsv",
    params:
        kmer_length=config["kmer_length"],
    conda:
        "../envs/kounta.yaml"
    threads: 8
    resources:
        mem_mb=4000,
        runtime=30
    shell:
        """
        ulimit -n 2048
        kounta --fofn {input} --kmer {params.kmer_length} --threads {threads} --ram 4 --out {output} 
        """

rule parse_kounta:
    input:
        "data/kmer_db/{gene}_kounta.tsv"
    output:
        "data/kmer_db/{gene}_kmer_query.fa",
        "data/kmer_db/{gene}_shared_kmer_list.txt",
        "data/kmer_db/{gene}_kmer_presence_absence.tsv"
    conda:
        "../envs/tidyverse.yaml"
    shell:
        """
        Rscript scripts/parse_kmer_counts.R {input} data/kmer_db/{wildcards.gene}
        """
