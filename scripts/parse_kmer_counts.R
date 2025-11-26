library(tidyverse)
library(glue)

# get command line arguments (kmer counts, output prefix)
args <- commandArgs(trailingOnly = TRUE)

# read in kounta output
t <- read_tsv(args[1])

# fix column names
t <- t %>% rename("kmer" = "#KMER") %>%  rename_with(~str_remove(., ".fasta"), everything())

# remove kmers present more than once in any allele
t <- t %>% filter(if_all(where(is.numeric), ~ . <=1))

# write updated kmer presence absence file
write_tsv(t, glue("{args[2]}_kmer_presence_absence.tsv"))

# identify kmers present in all alleles
conserved <- t %>% filter(if_all(where(is.numeric), ~ . == 1))
write_tsv(conserved %>% select(kmer), glue("{args[2]}_shared_kmer_list.txt"), col_names = F)

# write fasta file of query kmers for jellyfish
kmer_names <- t %>% select(kmer) %>% mutate(name = glue(">{kmer}"))
interleaved <- c(rbind(kmer_names$name, kmer_names$kmer))
file_connection <- file(glue("{args[2]}_kmer_query.fa"), open = "at")
writeLines(interleaved, file_connection)
close(file_connection)
