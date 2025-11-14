rule fsmlite_input:
    output:
        "data/kmer_db/{gene}_input_alleles.txt"
    params:
        allele_path=config["allele_path"]
    shell:
        """
        for f in {params.allele_path}/{wildcards.gene}/*.fasta; do id=$(basename "$f" .fasta); echo $id $f; done > {output}
        """

rule fsmlite:
    input:
        "data/kmer_db/{gene}_input_alleles.txt"
    output:
        kmers="data/kmer_db/{gene}_kmers.txt",
        tmp_meta=temp("{gene}_tmp.meta"),
        tmp_index=temp("{gene}_tmp.tmp")
    params:
        kmer_length=config["kmer_length"],
        total_alleles=103
    conda:
        "../envs/fsm-lite.yaml"
    shell:
        """
        mkdir -p kmers
        fsm-lite -l {input} -t {wildcards.gene}_tmp -m {params.kmer_length} -M {params.kmer_length} -s 1 -S {params.total_alleles} > {output.kmers}
        """

rule parse_fsmlite:
    input:
        "data/kmer_db/{gene}_kmers.txt"
    output:
        "data/kmer_db/{gene}_kmer_query.fa",
        "data/kmer_db/{gene}_shared_kmer_list.txt",
        "data/kmer_db/{gene}_unique_kmer_list.txt",
        "data/kmer_db/{gene}_kmer_presence_absence.tsv"
    params:
        total_alleles=103
    shell:
        """
        scripts/parse_allele_kmer_counts.py {input} {params.total_alleles} data/kmer_db/{wildcards.gene}
        """

