configfile: "config/config.yaml"

include: "rules/create_kmer_db.smk"

localrules: all

rule all:
    input:
        expand("data/kmer_db/{gene}_shared_kmer_list.txt", gene=["fimH", "sseL"]),
        expand("data/kmer_db/{gene}_unique_kmer_list.txt", gene=["fimH", "sseL"]),
        expand("data/kmer_db/{gene}_kmer_presence_absence.tsv", gene=["fimH", "sseL"])

