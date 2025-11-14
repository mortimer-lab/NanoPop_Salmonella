rule jellyfish_count:
    input:
        fastq=lambda wildcards: samples[wildcards.sample]["fastq"],
        query="data/kmer_db/{gene}_kmer_query.fa"
    output:
        "data/kmer_counts/{sample}_{gene}.jf"
    params:
        kmer_length=config["kmer_length"],
    conda:
        "../envs/jellyfish.yaml"
    threads: 8
    shell:
        """
        jellyfish count -m {params.kmer_length} -s 100M -C -t {threads} -o {output} --if {input.query} {input.fastq}
        """

rule jellyfish_dump:
    input:
        "data/kmer_counts/{sample}_{gene}.jf"
    output:
        "data/kmer_counts/{sample}_{gene}.txt"
    conda:
        "../envs/jellyfish.yaml"
    shell:
        """
        jellyfish dump -c {input} > {output}
        """
