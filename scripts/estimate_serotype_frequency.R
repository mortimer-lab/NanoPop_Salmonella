library(tidyverse)
library(stringi)
library(glue)

# get command line arguments (name of sample)
args <- commandArgs(trailingOnly = TRUE)

# read in kmer counts
fimH_counts <- read_tsv(glue("data/kmers_passing_threshold/{args[1]}_fimH.txt"), show_col_types = F)
sseL_counts <- read_tsv(glue("data/kmers_passing_threshold/{args[1]}_sseL.txt"), show_col_types = F)

# read in kmer presence absence matrices and add reverse complemented kmer
fimH_kmer_presence_absence <- read_tsv("data/kmer_db/fimH_kmer_presence_absence.tsv", show_col_types = F)
fimH_kmer_presence_absence <- fimH_kmer_presence_absence %>% mutate(reverse_complement = stri_reverse(chartr("ATCG", "TAGC", kmer)))
fimH_kmer_presence_absence <- fimH_kmer_presence_absence %>% relocate(reverse_complement, .after = kmer)


sseL_kmer_presence_absence <- read_tsv("data/kmer_db/sseL_kmer_presence_absence.tsv", show_col_types = F)
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
fimH_strains_present <- str_extract(fimH_alleles_present, "^[^-]+")

sseL_present <- sseL_strains %>% 
	map_lgl(check_strain_presence, counts = sseL_counts, presence_absence = sseL_kmer_presence_absence)
sseL_alleles_present <- sseL_strains[sseL_present]
sseL_strains_present <- str_extract(sseL_alleles_present, "^[^-]+")

strains_present <- sort(intersect(fimH_strains_present, sseL_strains_present))

# if Typhimurium and monophasic variant are both present, only include Typhimurium
if ("Typhimurium" %in% strains_present && "I4512i" %in% strains_present){
	strains_present <- strains_present[strains_present != "I4512i"]
}

# keep only the first allele from a "present" strain

fimH_alleles_present_filtered <- sort(str_subset(fimH_alleles_present, paste(paste0(strains_present, "-"), collapse='|')))
fimH_alleles_from_same_strain <- str_extract(fimH_alleles_present_filtered, "^[^-]+") %>% duplicated()
fimH_alleles_present_filtered <- fimH_alleles_present_filtered[!fimH_alleles_from_same_strain]

sseL_alleles_present_filtered <- sort(str_subset(sseL_alleles_present, paste(paste0(strains_present, "-"), collapse='|')))
sseL_alleles_from_same_strain <- str_extract(sseL_alleles_present_filtered, "^[^-]+") %>% duplicated()
sseL_alleles_present_filtered <- sseL_alleles_present_filtered[!sseL_alleles_from_same_strain]

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

fimH_unique_kmers <- get_unique_kmer(fimH_alleles_present_filtered, fimH_kmer_presence_absence)	
sseL_unique_kmers <- get_unique_kmer(sseL_alleles_present_filtered, sseL_kmer_presence_absence)

# get median counts for these unique kmers

get_median_counts <- function(kmer_list, kmer_counts){
	k_counts <- kmer_counts %>% filter(kmer %in% kmer_list) %>% pull(norm_count)
	return(median(k_counts))
}

fimH_median_counts <- fimH_unique_kmers %>% map_dbl(get_median_counts, kmer_counts = fimH_counts)
sseL_median_counts <- sseL_unique_kmers %>% map_dbl(get_median_counts, kmer_counts = sseL_counts)

# check for NAs in median counts; if they exist, remove counts for found serovars and rerun unique kmers/median counts

# define function to get minimum kmer count
get_minimum_counts <- function(kmer_list, kmer_counts){
	k_counts <- kmer_counts %>% filter(kmer %in% kmer_list) %>% pull(norm_count)
	return(min(k_counts))
}

# define function to update kmer counts with minimum count substracted from all present kmers

update_kmer_counts <- function(kmer_counts, presence_absence, min_strain, min_value){
	strain_kmers <- get_present_kmer(min_strain, presence_absence)
	new_kmer_counts <- kmer_counts %>% mutate(if_else(kmer %in% strain_kmers, norm_count - min_value, norm_count))
	return(new_kmer_counts)
}

# define functions to check if remaining alleles with unknown counts are in fact identical and unresolvable

check_identical_alleles <- function(unresolved_strain_list, presence_absence){
	filtered_patterns <- presence_absence %>% select(all_of(unresolved_strain_list))
	unique_patterns <- unique(as.list(filtered_patterns))
	# if there is only one unique presence absence pattern but multiple unresolved strains, return True
	return(length(unique_patterns) == 1 & length(unresolved_strain_list) > 1)
}

fimH_alleles_present_filtered_orig <- fimH_alleles_present_filtered
while (anyNA(fimH_median_counts)) {
	fimH_minimum_counts <- fimH_unique_kmers %>% map_dbl(get_minimum_counts, kmer_counts = fimH_counts)
	fimH_min_strain <- fimH_alleles_present_filtered[which.min(fimH_minimum_counts)]
	fimH_counts <- update_kmer_counts(fimH_counts, fimH_kmer_presence_absence, fimH_min_strain, min(fimH_minimum_counts))
	fimH_alleles_present_filtered <- fimH_alleles_present_filtered[fimH_alleles_present_filtered != fimH_min_strain] 
	fimH_unique_kmers <- get_unique_kmer(fimH_alleles_present_filtered, fimH_kmer_presence_absence)	
	fimH_median_counts_new <- fimH_unique_kmers %>% map_dbl(get_median_counts, kmer_counts = fimH_counts)
	# update counts with new values
	for (a in fimH_alleles_present_filtered){
		orig_index <- match(a, fimH_alleles_present_filtered_orig)
		new_index <- match(a, fimH_alleles_present_filtered)
		fimH_median_counts[orig_index] <- fimH_median_counts_new[new_index]
	}
	if (check_identical_alleles(fimH_alleles_present_filtered_orig[is.na(fimH_median_counts)], fimH_kmer_presence_absence)){
		break
	}
}

sseL_alleles_present_filtered_orig <- sseL_alleles_present_filtered
while (anyNA(sseL_median_counts)) {
	sseL_minimum_counts <- sseL_unique_kmers %>% map_dbl(get_minimum_counts, kmer_counts = sseL_counts)
	sseL_min_strain <- sseL_alleles_present_filtered[which.min(sseL_minimum_counts)]
	sseL_counts <- update_kmer_counts(sseL_counts, sseL_kmer_presence_absence, sseL_min_strain, min(sseL_minimum_counts))
	sseL_alleles_present_filtered <- sseL_alleles_present_filtered[sseL_alleles_present_filtered != sseL_min_strain] 
	sseL_unique_kmers <- get_unique_kmer(sseL_alleles_present_filtered, sseL_kmer_presence_absence)
	sseL_median_counts_new <- sseL_unique_kmers %>% map_dbl(get_median_counts, kmer_counts = sseL_counts)
	# update counts with new values
	for (a in sseL_alleles_present_filtered){
		orig_index <- match(a, sseL_alleles_present_filtered_orig)
		new_index <- match(a, sseL_alleles_present_filtered)
		sseL_median_counts[orig_index] <- sseL_median_counts_new[new_index]
	}
	if (check_identical_alleles(sseL_alleles_present_filtered_orig[is.na(sseL_median_counts)], sseL_kmer_presence_absence)){
		break
	}
}


# read in median counts from conserved kmers

fimH_conserved_counts <- scan(glue("data/conserved_kmer_count/{args[1]}_fimH.txt"), what = numeric())
sseL_conserved_counts <- scan(glue("data/conserved_kmer_count/{args[1]}_sseL.txt"), what = numeric())

# make output data frame

median_counts <- data.frame(strain = strains_present, fimH_counts = fimH_median_counts, sseL_counts = sseL_median_counts)

median_counts <- median_counts %>% 
	mutate(fimH_known_proportion = fimH_counts/sum(fimH_counts)) %>%
	mutate(sseL_known_proportion = sseL_counts/sum(sseL_counts)) %>%
	mutate(average_known_proportion = (fimH_known_proportion+sseL_known_proportion)/2)

# estimate proportion of reads covered by known strains

fimH_unknown <- 100*pmax(0, 1 - sum(median_counts$fimH_counts)/fimH_conserved_counts)
sseL_unknown <- 100*pmax(0, 1 - sum(median_counts$sseL_counts)/sseL_conserved_counts)

if (anyNA(c(fimH_unknown, sseL_unknown))){
	unknown_sample_description <- glue("{args[1]}\tsample contains two strains with identical alleles in fimH or sseL : unknowns cannot be identified")
} else if (fimH_unknown < 0.985 && sseL_unknown < 0.985) {
	unknown_sample_description <- glue("{args[1]}\tall strains represented in database")
} else {
	unknown_sample_description <- glue("{args[1]}\tpotential unknown strain: {fimH_unknown}% of fimH reads do not match known allele, {sseL_unknown}% of sseL reads do not match known allele")
}

file_connection <- file("data/strain_proportion_estimates/unknown_strain_estimates.txt", open = "a")
writeLines(unknown_sample_description, file_connection) 
close(file_connection)

write_tsv(median_counts, glue("data/strain_proportion_estimates/{args[1]}.txt"))
