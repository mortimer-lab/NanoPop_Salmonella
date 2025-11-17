import pandas as pd

configfile: "config/config.yaml"

include: "rules/create_kmer_db.smk"
include: "rules/count_kmers.smk"
include: "rules/normalize_kmers.smk"

localrules: all,normalize_neg

samples = pd.read_table(config["samples"]).set_index("sample").to_dict(orient="index")
samples_without_negative = list(samples.keys())
samples_without_negative.remove(config["negative_control_sample"])

rule all:
    input:
        expand("data/negative_control_normalized/{sample}_{gene}.txt", gene=["fimH", "sseL"], sample=samples_without_negative)

