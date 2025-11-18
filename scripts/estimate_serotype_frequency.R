library(tidyverse)
library(stringi)
library(glue)

# get command line arguments (name of sample)
args <- commandArgs(trailingOnly = TRUE)

# read in kmer counts
fimH_counts <- read_tsv(glue("data/kmers_passing_threshold/{args[1]}_fimH.txt"))
sseL_counts <- read_tsv(glue("data/kmers_passing_threshold/{args[1]}_sseL.txt"))

# read in kmer presence absence matrices and add reverse complemented kmer
fimH_kmer_presence_absence <- read_tsv("data/kmer_db/fimH_kmer_presence_absence.tsv")
fimH_kmer_presence_absence <- fimH_kmer_presence_absence %>% rename("kmer" = "...1")
fimH_kmer_presence_absence <- fimH_kmer_presence_absence %>% mutate(reverse_complement = stri_reverse(chartr("ATCG", "TAGC", kmer)))
fimH_kmer_presence_absence <- fimH_kmer_presence_absence %>% relocate(reverse_complement, .after = kmer)


sseL_kmer_presence_absence <- read_tsv("data/kmer_db/sseL_kmer_presence_absence.tsv")
sseL_kmer_presence_absence <- sseL_kmer_presence_absence %>% rename("kmer" = "...1")
sseL_kmer_presence_absence <- sseL_kmer_presence_absence %>% mutate(reverse_complement = stri_reverse(chartr("ATCG", "TAGC", kmer)))
sseL_kmer_presence_absence <- sseL_kmer_presence_absence %>% relocate(reverse_complement, .after = kmer)

# get strain names 
fimH_strains <- tail(colnames(fimH_kmer_presence_absence), -2)
sseL_strains <- tail(colnames(sseL_kmer_presence_absence), -2)


# function to detect presence of strains

check_strain_presence <- function(strain_name, counts, presence_absence){
	k <- presence_absence %>% 
		filter(.data[[strain_name]] == 1) %>% 
		select(kmer, reverse_complement) %>% 
		left_join(counts, by = join_by("kmer" == "kmer")) %>% 
		left_join(counts, by = join_by("reverse_complement" == "kmer")) %>% 
		mutate(across(everything(), ~replace_na(.,0))) %>%
		mutate(final_count = pmax(norm_count.x, norm_count.y, na.rm = T))
	minimum_serotype_kmer_count <- min(k$final_count)
	return(minimum_serotype_kmer_count > 0)
}

# check for presence of strains via fimH and sseL
fimH_present <- fimH_strains %>% 
	map_lgl(check_strain_presence, counts = fimH_counts, presence_absence = fimH_kmer_presence_absence)

fimH_alleles_present <- fimH_strains[fimH_present]
sseL_present <- sseL_strains %>% 
	map_lgl(check_strain_presence, counts = sseL_counts, presence_absence = sseL_kmer_presence_absence)

sseL_alleles_present <- sseL_strains[sseL_present]

strains_present <- intersect(fimH_alleles_present, sseL_alleles_present)

# get unique fimH and sseL kmers for each strain
get_present_kmer <- function(strain, presence_absence){
	kmers <- presence_absence %>% 
		filter(.data[[strain]] == 1) %>% 
		pull(kmer)
	reverse_comp <- presence_absence %>% 
		filter(.data[[strain]] == 1) %>% 
		pull(reverse_complement)
	all <- c(kmers, reverse_comp)
}

get_unique_kmer <- function(strain_list, presence_absence){
	strain_kmers <- strain_list %>% map(get_present_kmer, presence_absence = presence_absence)
	unique_kmers <- lapply(1:length(strain_kmers), function(n) setdiff(strain_kmers[[n]], unlist(strain_kmers[-n])))
	return(unique_kmers)
}

fimH_unique_kmers <- get_unique_kmer(strains_present, fimH_kmer_presence_absence)	
sseL_unique_kmers <- get_unique_kmer(strains_present, sseL_kmer_presence_absence)

# get median counts for these unique kmers

get_median_counts <- function(kmer_list, kmer_counts){
	k_counts <- kmer_counts %>% filter(kmer %in% kmer_list) %>% pull(norm_count)
	return(median(k_counts))
}

fimH_median_counts <- fimH_unique_kmers %>% map_dbl(get_median_counts, kmer_counts = fimH_counts)
sseL_median_counts <- sseL_unique_kmers %>% map_dbl(get_median_counts, kmer_counts = sseL_counts)

median_counts <- data.frame(strain = strains_present, fimH_counts = fimH_median_counts, sseL_counts = sseL_median_counts)

median_counts <- median_counts %>% 
	mutate(fimH_proportion = fimH_counts/sum(fimH_counts)) %>%
	mutate(sseL_proportion = sseL_counts/sum(sseL_counts)) %>%
	mutate(average_proportion = (fimH_proportion+sseL_proportion)/2)

write_tsv(median_counts, glue("data/strain_proportion_estimates/{args[1]}.txt"))
