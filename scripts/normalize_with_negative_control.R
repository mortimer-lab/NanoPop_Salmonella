library(tidyverse)
library(glue)

# get command line arguments (names of gene, sample, negative control)
args <- commandArgs(trailingOnly = TRUE)

# read in kmer counts from sample and negative control
counts <- read_delim(glue("data/kmer_counts/{args[2]}_{args[1]}.txt"), col_names = c("kmer", "count"))
neg <- read_delim(glue("data/kmer_counts/{args[3]}_{args[1]}.txt"), col_names = c("kmer", "negative_count"))

# join dataframes
joined <- counts %>% left_join(neg)

# remove counts found in negative control, if count in negative control is greater than count in sample, set count to zero
joined <- joined %>% mutate(norm_count = if_else(count > negative_count, count-negative_count, 0))

# save output to TSV
write_tsv(joined %>% select(kmer, norm_count), glue("data/negative_control_normalized/{args[2]}_{args[1]}.txt"))
