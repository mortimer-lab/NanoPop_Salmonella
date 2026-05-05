# NanoPop_Salmonella
Snakemake pipeline for automated estimation of Salmonella serovar frequencies within a mixed population using analysis of amplicon sequencing of fimH and sseL with Nanopore

## Pipeline Requirements
This pipeline requires conda for software installation and Snakemake. If you are running this pipeline on Sapelo2 or another cluster with slurm, you also need the snakemake executor plugin for slurm.

### Installing snakemake
1. Install conda or mamba. We recommend [miniforge](https://conda-forge.org/download/).
2. Install [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

#### Required only if you are using slurm:

3. Install [Snakemake executor plugin: slurm](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) in your snakemake conda environment.

## Pipeline Components

### envs
YAML files describing conda environments required for this pipeline

### config
Configuration files specific to your analysis. Template files have been included in this directory. The components are described below:

#### config.yaml
```
samples: "config/samples.tsv"                 # path to samples file

allele_path: "alleles"                        # path to alleles database

kmer_length: 15                               # kmer length to create kmer presence/absence database

negative_control_sample: "negative"           # name of negative control sample

kmer_frequency_threshold: 0.002               # minimum frequency of strain specific kmers compared to conserved kmers to consider a strain present
```
#### samples.tsv
The sample description file should have sample names in the first column and paths to fastqs in the second column. The header row is required. The sample name for the negative control should match what is listed in the `config.yaml` file.
```
sample	fastq
sample1	data/fastqs/sample1.fastq
sample2	data/fastqs/sample2.fastq
sample3	data/fastqs/sample3.fastq
negative	data/fastqs/negative.fastq
```

### rules
Snakemake files describing rules for pipeline components. These should not need to be edited by the user. However, if you find that a particular rule does not have sufficient resources for your analyses. The resources block can be added or edited. For example:

```
rule kounta:
    input:
        "data/kmer_db/{gene}_input_alleles.txt"
    output:
        kmers="data/kmer_db/{gene}_kounta.tsv",
    params:
        kmer_length=config["kmer_length"],
    conda:
        "../envs/kounta.yaml"
    threads: 8
    resources:
        mem_mb=4000,                                       # expected memory needed (in MB)
        runtime=30                                         # expected time needed (in minutes)
    shell:
        """
        ulimit -n 2048
        kounta --fofn {input} --kmer {params.kmer_length} --threads {threads} --ram 4 --out {output} 
        """
```

### scripts
Custom R scripts for data analysis. The use of these scripts is automated by the Snakemake pipeline. Their individual purpose and arguments are described below. 

#### apply_threshold.R
This script applies the specified kmer threshold to remove very low frequency kmers.

Usage: `Rscript apply_threshold.R gene_name sample_name threshold`

#### estimate_serotype_frequency.R
This script estimates the frequency of different serotypes using normalized kmer counts from *fimH* and *sseL*.

Usage: `Rscript estimate_serotype_frequency.R sample_name`

#### normalize_with_negative_control.R
This script removes kmer counts found in the negative control sample.

Usage: `Rscript normalize_with_negative_control.R gene_name sample_name negative_control_name`

#### parse_kmer_counts.R
This script reads the kmer counts file output by kounta and identifies conserved kmers (present in all alleles), creates a cleaned kmer presence absence file, and writes a fasta file of query kmers for jellyfish.

Usage: `Rscript parse_kmer_counts.R kmer_counts_file output_prefix`

### Snakefile
Main Snakefile that imports rules and defines desired output. This file should not need to be edited by the user.

## Required Input

### alleles
A database of alleles for each gene is required. These alleles should be in a directory with the name of the gene, and each allele should be in a separate fasta file. It is required that the allele be named with the strain/serotype name followed by a dash and then a specific genome designation from which the allele came from. For example:

```
alleles/

fimH/
Agona-GCA_006890265.1.fasta        KentuckyI-GCA_008516425.1.fasta      Poona-GCA_008480745.1.fasta
Anatum-GCA_006839465.1.fasta       KentuckyI-GCA_008638905.1.fasta      Potsdam-GCA_006304605.1.fasta

sseL/
Agona-GCA_006890265.1.fasta        KentuckyI-GCA_008516425.1.fasta      Poona-GCA_008480745.1.fasta
Anatum-GCA_006839465.1.fasta       KentuckyI-GCA_008638905.1.fasta      Potsdam-GCA_006304605.1.fasta
```

The directory containing the alleles can be specified in the `config/config.yaml' file. An allele database has been included with this repository, but additional alleles can be added.

### fastqs for samples and negative control
fastq files with sequencing output for each sample and a negative control is required. The path to these fastqs should be provided in a tab-separated file with sample names in the first column and fastq paths in the second column (see example in `config/samples.tsv`). The name of your sample file can be specificied in `config/config.yaml`. Additionally, the sample name for your negative control should also be listed in this file.

## Running the pipeline

1. Clone this repository.
2. Update `config/config.yaml` and `config/samples.tsv` with information specific to your dataset.
3. Make sure all desired output is listed under `rule all:` in Snakefile
4. Run snakemake using `snakemake`. If you would like to run on a cluster using slurm, specify the profile: `snakemake --workflow-profile profiles/slurm`.

## Pipeline Output
All pipeline intermediate files and output is saved in the `data` directory.

### conserved_kmer_count
Counts of conserved kmers in each sample for *fimH* and *sseL*.

### kmer_counts
Counts of query kmers in each sample for *fimH* and *sseL*.

### kmer_db
Kmer presence absence database for *fimH* and *sseL*. Includes both a presence absence matrix and a query fasta for jellyfish.

### kmers_passing_threshold
Counts of query kmers that pass specified frequency thresholds in each sample for *fimH* and *sseL*.

### negative_countrol_normalized
Counts of query kmers after normalization based on negative control counts in each sample for *fimH* and *sseL*.

### strain_proportion_estimates
Tab-separated text files for each sample showing serotype-specific *fimH* and *sseL* kmer counts, estimated proportions of known serotypes based on *fimH* and *sseL* kmer counts separately, and estimated proportions of known serotypes based on the average of *fimH*- and *sseL*-specific estimates.

There is also a file called *unknown_strain_estimates.txt* that indicates if there is evidence of an unknown serotype in the sample.

## Final Output File Descriptions
Final output files can be found in `data/strain_proportion_estimates/`.

### Example sample-specific output
In cases where serotypes have kmer presence/absence patterns that cannot be distinguished for either *fimH* or *sseL*, counts will be reported as "NA" for those serotypes and a known proportion will not be estimated.

```
strain  fimH_counts     sseL_counts     fimH_known_proportion   sseL_known_proportion   average_known_proportion
Infantis        25899   11196   0.9380636748886233      0.9550456367823936      0.9465546558355085
KentuckyI       1710    527     0.06193632511137673     0.044954363217606416    0.05344534416449157
```
### Example unknown strain estimate output
```
sample1   all strains represented in database
sample2   all strains represented in database
sample3   potential unknown strain: 2.9% of fimH reads do not match known allele, 4.7% of sseL reads do not match known allel
```



