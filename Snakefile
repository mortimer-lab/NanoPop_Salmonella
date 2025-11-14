import pandas as pd

configfile: "config/config.yaml"

include: "rules/create_kmer_db.smk"
include: "rules/count_kmers.smk"

localrules: all

samples = pd.read_table(config["samples"]).set_index("sample").to_dict(orient="index")

rule all:
    input:
        expand("data/kmer_counts/{sample}_{gene}.txt", gene=["fimH", "sseL"], sample=samples.keys())

