# ERMA - epicPCR Resistome-Microbiome Analyzer

[![Snakemake CI](https://github.com/IKIM-Essen/ERMA/actions/workflows/snakemake-ci.yml/badge.svg)](https://github.com/IKIM-Essen/ERMA/actions/workflows/snakemake-ci.yml)

This pipeline is designed to process sequencing reads obtained from epicPCR experiments, linking antimicrobial resistance (AMR) genes with 16S rRNA genes. The pipeline uses **Snakemake** to manage the workflow and integrates tools for downloading necessary databases, running sequence alignments, filtering, and generating visual reports.

## Features
- Downloads and prepares SILVA and CARD databases for use in the analysis.
- Converts raw FASTQ sequencing files to FASTA format.
- Performs BLAST sequence alignments against both CARD (AMR genes) and SILVA (16S rRNA).
- Integrates BLAST results to identify linked AMR and microbial markers.
- Generates filtered results based on alignment quality and similarity thresholds.
- Produces graphical outputs such as genus abundance plots, boxplots, and tables.
- Generates an HTML report summarizing the analysis.

## Prerequisites

### Install Snakemake
The pipeline requires **Snakemake** with support for conda environments. You can install it using conda (via Miniconda or Anaconda):

```bash
mamba install -c conda-forge -c bioconda snakemake --yes; mamba install -c bioconda snakemake-storage-plugin-fs --yes
```

Alternatively, follow the official [Snakemake installation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) guide for more options.

Install Dependencies
This pipeline uses conda environments to manage its dependencies. Snakemake will automatically create and manage these environments when run with the `--use-conda` flag (default in profile).

## Usage Instructions

Clone the repository: First, clone the pipeline repository to your local machine:

```bash
git clone https://github.com/your-username/ERMA.git
cd ERMA
```
Prepare Data Folder: You need to place your raw sequencing files (fastq.gz format) in the data/fastq/ directory. This folder must exist before running the pipeline.

Modify the Config File: Open the config/config.yaml file and change the base_dir parameter to the base directory where the pipeline is located. The config file should look like this:

```yaml
runname: "ERMA_runname123"
base_dir: "/path/to/your/ERMA"

min_similarity: "0.8" # threshold to filter blast hits by percentage identity

silva:
  download-path-seq: "path/to/silva_db"
  download-path-tax: "path/to/silva_taxmap"

card:
  download-path: "https://card.mcmaster.ca/download/0/broadstreet-v3.3.0.tar.bz2"

num_parts: 1 # number of subfiles the fastqs are split into
max_threads: 16 # 

similarity_search_mode: "extensive" # Put here fast or extensive for search on only forward or both strands

seq_tech: "ONT"
```

Replace /path/to/your/project with the actual path to your local pipeline directory.

### Illumina input

Skip this when analyzing single-end or already merged paired-end reads.

When analyzing paired-end reads, they must be merged before starting the pipeline. For this ERMA provides a module. Load your paired end FASTQ-files in the data/fastq folder and execute from the folder where also the Snakefile is

1. Copy the paired-end fastq.gz files in data/fastq
2. Set the seq_tech parameter to "Illumina"
3. Run:

```bash
snakemake prepare-fastqs
```

### ONT input

Skip this when the input files are already in fastq.gz format.

When starting with raw ONT output, this command can be used to demultiplex the samples according to the respective barcodes. For that, following must be prepared:

1. Enter the desired Sample name - barcode combinations in the file "barcode-rename.csv" which can be found at data/ONT
2. Copy the fastq_pass folder of the targeted ONT run in data/ONT or change the respective path in the config file
3. Define the target fragment lengths parameter in the config file - everything outside the filter intervall parameter gets filtered out beforehand
4. Set the seq_tech parameter to "ONT"
5. Run:

```bash
snakemake prepare-fastqs
```

### Run Pipeline

To run the pipeline in the default mode, a profile with snakemake specific parameters is provided. Therefore, the only neccessary command to start is:

```bash
snakemake --profile profile
```

## Additional Notes

The pipeline is designed to handle large sequencing datasets in parallel, so it's recommended to run it on a machine with sufficient computational resources. However, to run the pipeline on machines with less resources, it is recommended to split the fastq files or the tables in smaller chunks to prevent the RAM to overflow. This can be done by increase num_parts in the config file. However, the higher the number of parts per FASTQ file, the higher the chance some blast results for the same read will be split and lost.
If any errors occur during the pipeline run, Snakemake will provide detailed logs, allowing you to debug and troubleshoot any issues. You are most welcome to create an Issue when running into problems.

License

This project is licensed under the MIT License.
