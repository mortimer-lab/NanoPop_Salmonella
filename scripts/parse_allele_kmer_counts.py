#!/usr/bin/env python

import argparse
import pandas as pd
from collections import defaultdict
from collections import Counter

def get_args():
    parser = argparse.ArgumentParser(description='Parse fsm-lite output')
    parser.add_argument("kmers", help="fsm-lite output file")
    parser.add_argument("total_sequences", help="Total alleles", type=int)
    parser.add_argument("output_prefix", help="prefix for output files")
    return parser.parse_args()

def parse_file(kmers, total, prefix):
    kmer_dict = defaultdict(list)
    with open(kmers, "r") as infile:
        for line in infile:
            line = line.strip().split(" | ")
            kmer = line[0]
            alleles = []
            for a in line[1].split():
                count = a.split(":")[1]
                if count != "1": # only include single copy kmers
                    break
                else:
                    name = a.split(":")[0]
                    kmer_dict[kmer].append(name)
    with open(f"{prefix}_shared_kmer_list.txt", "w") as shared_file:
        with open(f"{prefix}_unique_kmer_list.txt", "w") as unique_file:
            with open(f"{prefix}_kmer_query.fa", "w") as jellyfish_query:
                for k in kmer_dict:
                    jellyfish_query.write(f">{k}\n{k}\n")
                    if len(kmer_dict[k]) == total:
                        shared_file.write(k + "\n")
                    else:
                        unique_file.write(k + "\n")
    df = pd.DataFrame({k:Counter(v) for k,v in kmer_dict.items()}).T.fillna(0).astype(int)
    df.to_csv(f"{prefix}_kmer_presence_absence.tsv", sep="\t")

def main():
    args = get_args()
    parse_file(args.kmers, args.total_sequences, args.output_prefix)

if __name__ == "__main__":
    main()
