library(tidyverse)
library(stringi)
library(glue)

# get command line arguments (names of gene, sample, threshold proportion)
args <- commandArgs(trailingOnly = TRUE)

# read in kmer counts and list of conserved kmers
kmer_counts <- read_tsv(glue("data/negative_control_normalized/{args[2]}_{args[1]}.txt"))
conserved_kmers <- read_tsv(glue("data/kmer_db/{args[1]}_shared_kmer_list.txt"), col_names = c("kmer"))

# reverse complement conserved kmers
conserved_kmers <- conserved_kmers %>% mutate(reverse_complement = stri_reverse(chartr("ATCG","TAGC", kmer))
)

# get median of conserved kmers and apply threshold
kmer_counts_conserved <- kmer_counts %>% filter(kmer %in% conserved_kmers$kmer | kmer %in% conserved_kmers$reverse_complement)
median_conserved_count <- median(kmer_counts_conserved$norm_count)
threshold <- median_conserved_count*as.numeric(args[3])

# filter all kmer counts to those greater than threshold
kmer_counts_filtered <- kmer_counts %>% filter(norm_count > threshold)
write_tsv(kmer_counts_filtered, glue("data/kmers_passing_threshold/{args[2]}_{args[1]}.txt"))

