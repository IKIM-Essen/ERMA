# ERMA - epicPCR Resistome-Microbiome Analyzer

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Snakemake CI](https://github.com/IKIM-Essen/ERMA/actions/workflows/snakemake-ci.yml/badge.svg)](https://github.com/IKIM-Essen/ERMA/actions/workflows/snakemake-ci.yml)

This **Snakemake**-based pipeline processes sequencing reads from epicPCR experiments to link antimicrobial resistance (AMR) genes with genus-specific 16S rRNA gene sequences from prokaryotic cells. It orchestrates all workflow steps, including database downloads, sequence alignment, read filtering, and the generation of visual reports, ensuring reproducibility and streamlined analysis.

## Features
- Downloads and prepares SILVA and CARD databases for use in the analysis.
- Converts raw FASTQ sequencing files to FASTA format.
- Performs Diamond/Usearch sequence alignments against both CARD (AMR genes) and SILVA (16S rRNA).
- Integrates similarity search results to identify linked AMR and microbial markers.
- Generates filtered results based on alignment quality and similarity thresholds.
- Produces graphical outputs such as genus abundance plots, boxplots, and tables.
- Generates an HTML report summarizing the analysis.

## Prerequisites

### Minimum Installation

This pipeline is possible to deploy with snakedeploy. We recommend to only use this with already present experience with snakemake and snakedeploy. If this is not the case, we recommend to use the full installation guide.

For usage:
1. Install snakedeploy and snakemake≥v9 (if necessary)
2. deploy minimum workflow distribution (with desired destination)
```bash
snakedeploy deploy-workflow https://github.com/IKIM-Essen/ERMA --dest . --branch snakedeploy
```
3. prepare setup
```bash
mkdir -p data/fastq
cp "files-to-analyze" data/fastq
```
4. start pipeline (change config if necessary)
```bash
snakemake --cores all --sdm conda
```

### Full Installation
The pipeline requires **Snakemake** (≥ Version 9.0) with support for conda environments. You can install it using conda/mamba (via Miniconda or Anaconda):

```bash
mamba create -n snakemake9 -c conda-forge -c bioconda snakemake=9 --yes
mamba activate snakemake9
```

Alternatively, follow the official [Snakemake installation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) guide for more documentation and options.

Install Dependencies: This pipeline uses conda environments to manage its dependencies. Snakemake will automatically create and manage these environments when run with the `--use-conda` flag (default in profile).

## Usage Instructions

First, clone the pipeline repository to your local machine:

```bash
git clone https://github.com/IKIM-essen/ERMA.git
cd ERMA
```
Prepare Data Folder: You need to place your raw sequencing files (fastq.gz format) in the data/fastq/ directory or change this path to the desired directory.

Modify the Config File: Open the config/config.yaml file and change the base_dir parameter to the base directory where the pipeline is located. The config file should look like this:

```yaml
runname: "ERMA_runname123"
# setting up base directory and location of input and output. Generally, no changes needed here.
base_dir: "."
fastq_dir: "data/fastq" # copy target fastq.gz files in ERMA/data/fastq or change this path
outdir: "results" # Output directory of the final report

min_similarity: "0.8" # threshold to filter blast hits by percentage identity
min_abundance: "0.01" # genera with lower abundance will be binned as "Other" in stacked bar abundance plot

silva:
  download_path_seq: "path/to/silva_db"
  download_path_tax: "path/to/silva_taxmap"

card:
  download_path: "https://card.mcmaster.ca/download/0/broadstreet-v3.3.0.tar.bz2"

num_parts: 1 # number of subfiles the fastqs are split into
max_threads: 16 

similarity_search_mode: "full" # Put here "test" or "full" for strand/s to be included in the similarity search

### Preprocessing ###
# if data is already in format 'one fastq.gz per sample', this section can be ignored

seq_tech: "Illumina" # Put here "Illumina" or "ONT" befor using rule prepare_fastqs

# In case of Demultiplexing ONT data, provide information for this section
ONT:
  fastq_pass_path: "data/ONT/fastq_pass" # copy your fastq_pass folder here
  sample_name_path: "data/ONT/barcode-rename.csv" # change this file with your barcode-sample name combinations
  target_fragment_length: 1250 # Length of the theoretical fragment after nested PCR
  filter_intervall: 0.1 # +/- Intervall used to filter too large/small fragments; 0.1 filters in a +/- 10% intervall
```

### Illumina input

Skip this when analyzing single-end or already merged paired-end reads.

When analyzing paired-end reads, they must be merged before starting the pipeline. For this, ERMA provides a module. Follow the steps:

1. Copy the to-be-merged paired-end fastq.gz files in data/fastq
2. Set the seq_tech parameter to "Illumina" in config/config.yaml (if not already)
3. Run from the ERMA root folder:

```bash
snakemake prepare_fastqs
```

This will execute a bash script with the essential step of merging R1 and R2 read with NGmerge

### ONT input

Skip this when the input files are already in "one fastq.gz file per sample" - format.

When starting with raw ONT output, this routine can be used to demultiplex the samples according to the respective barcodes. For that, following must be prepared:

1. Enter the desired sample name - barcode combinations in the file "barcode-rename.csv" which can be found at data/ONT
2. Copy the fastq_pass folder of the targeted ONT run in data/ONT or change the respective path in the config file
3. Define the target fragment lengths parameter in the config file - everything outside the filter intervall parameter gets filtered out beforehand
4. Set the seq_tech parameter to "ONT" in config/config.yaml (if not already)
5. Run from the ERMA root folder:

```bash
snakemake prepare_fastqs
```

This will execute a python script that filters and demultiplex the ONT output in one fastq file per sample.

### Sample naming

To ensure a robust running of the pipeline, sample names are restricted to only contain:
- Uppercase and/or lowercase letters
- digits 
- underscores

Sample names containing other special characters ("-" "," "." "/" "$" "@" etc.) will result in an error

### Run Pipeline

To run the pipeline in the default mode, a profile with snakemake specific parameters is provided. Therefore, the only neccessary command to start is:

```bash
snakemake --profile profile
```

### Testing

For testing the workflow you can copy the provided dummy data:

```
cp .github/data/fastq/test_epic_data.fastq.gz data/fastq/
```

In this case, the similarity search mode in the config file can be changed to "test"

## Additional Notes

The pipeline is designed to handle large sequencing datasets in parallel, so it's recommended to run it on a machine with sufficient computational resources. However, to run the pipeline on machines with less resources, it is recommended to split the fastq files or the tables in smaller chunks to prevent RAM  overflow. This can be done by increase num_parts in the config file. However, the higher the number of parts per FASTQ file, the higher the chance some blast results for the same read will be split and lost.
If any errors occur during the pipeline run, Snakemake will provide detailed logs, allowing you to debug and troubleshoot any issues. You are most welcome to create an Issue when running into problems.

License

This project is licensed under the MIT License.
