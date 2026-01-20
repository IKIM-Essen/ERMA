# ERMA - epicPCR Resistome-Microbiome Analyzer

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Snakemake CI](https://github.com/IKIM-Essen/ERMA/actions/workflows/snakemake-ci.yml/badge.svg)](https://github.com/IKIM-Essen/ERMA/actions/workflows/snakemake-ci.yml)

This **Snakemake**-based pipeline processes sequencing reads from epicPCR experiments to link antimicrobial resistance (AMR) or other target genes with genus-specific 16S rRNA gene sequences from prokaryotic cells. It orchestrates all workflow steps, including database downloads, sequence alignment, read filtering, and the generation of visual reports, ensuring reproducibility and streamlined analysis.

## Feature overview
- Auto-downloads and prepares SILVA and CARD databases for use in the analysis
- Auto-downloads user-specified targets (UniRef )and adds it to the CARD database
- Performs Diamond/Usearch sequence alignments against both databases
- Integrates similarity search results to identify linked AMR and microbial markers
- Generates filtered results based on alignment quality and similarity thresholds
- Produces graphical outputs such as genus abundance plots, boxplots, and tables
- Generates an HTML report summarizing the analysis

## Prerequisites

### Minimum Installation

This pipeline is possible to deploy with snakedeploy. We recommend to only use this with already present experience with snakemake and snakedeploy. If this is not the case, we recommend to use the full installation guide.

For usage:
1. Make sure snakedeploy and snakemake≥v9 is installed
2. deploy minimum workflow distribution (with desired destination)
```bash
snakedeploy deploy-workflow https://github.com/IKIM-Essen/ERMA --dest . --branch snakedeploy
```
3. prepare setup
```bash
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

### Preparing Databases

ERMA relies on two primary reference databases: SILVA for taxonomic assignment and CARD for antimicrobial resistance gene detection. Both databases are handled directly through the workflow configuration and can be fetched automatically or provided locally, depending on user requirements.

#### Automatic Download (Default)

By default, ERMA downloads both databases using the URLs defined in the configuration file:

```yaml
download_silva: "https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
download_card: "https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2"
```

When these entries remain set to URLs, the workflow retrieves and processes the databases automatically during execution. No additional preparation steps are required. Once downloaded, the workflow uses the present files for reiterations.

#### Using Local Database Files

If the user prefers handling the downloads manually or already has the required databases available on disk, the workflow can operate entirely from local paths. To do so, replace the URLs in the configuration file with absolute file paths:

```yaml
download_silva: "/path/to/local/SILVA_138.2.fasta.gz"
download_card: "/path/to/local/broadstreet-v4.0.1.tar.bz2"
```

Once local paths are set, ERMA will skip remote downloads and use the supplied files without modification of the pipeline structure.

#### Adding UniRef Targets of Interest (Optional)

ERMA provides an optional mechanism to incorporate UniRef clusters representing additional genes or functions relevant to the analysis. This process is fully automated and controlled through the following section of the configuration file:

```yaml
add_uniref_targets: False
uniprot_cluster: "100"
uniprot_targets: ["int1", "inti1", "class_1_integron"]
max_entry_count: 1000
low_freq_threshold: 0.01
```

To activate UniRef integration, set:

```yaml
add_uniref_targets: True
```

Then adjust the uniprot_targets list to include any desired search terms or gene names. ERMA downloads the corresponding UniRef entries, filters them according to the configured uniprot_cluster level (e.g., UniRef100), and integrates them directly into the CARD-derived reference database used for similarity searches.
This option enables users to extend the AMR and functional screening capabilities of ERMA without requiring manual reference curation.

#### Manual Addition of New Targets (Alternative Workflow)

Users who prefer to add custom targets manually can follow this approach:
1. Retrieve the protein FASTA sequences for the target gene(s) from UniProt, CARD, or any external source
2. Append the new sequences to the CARD protein FASTA file used by ERMA (protein_fasta_protein_homolog_model.fasta)
3. Rezip the beforehand unzipped files to the same format as before
4. Modify the workflow configuration to point to the updated local database files
This approach offers full control over the source, format, and curation of custom targets while keeping the surrounding pipeline unchanged

## Pipeline Usage

First, clone the pipeline repository to your local machine:

```bash
git clone https://github.com/IKIM-essen/ERMA.git
cd ERMA
```
Prepare Data Folder: You need to place your raw sequencing files (fastq.gz format) in the data/fastq/ directory or change this path to the desired directory.

Modify the Config File: Open the config/config.yaml file and change the base_dir parameter to the base directory where the pipeline is located. The most relevant parameters for the standard user are:

```yaml
runname: "runname123"
min_similarity: "0.9" # threshold to filter blast hits by percentage identity
num_parts: 1 # number of subfiles the fastqs are split into (prevents crashing with very large input)
max_threads: 16 
seq_tech: "Illumina" # Put here "Illumina" or "ONT" according to used technology (only important for using rule prepare_fastqs)
```

### Illumina input

Skip this when analyzing single-end or already merged paired-end reads.

When analyzing paired-end reads, they must be merged before starting the pipeline. For this, ERMA provides a module. Follow the steps:

1. Copy the to-be-merged paired-end fastq.gz files in data/fastq
2. Set the seq_tech parameter to "Illumina" in config/config.yaml (if not already)
3. Run from the ERMA root folder:

```bash
snakemake prepare_fastqs --cores N
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
snakemake prepare_fastqs --cores N
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

For testing the workflow the user can copy the provided dummy data:

```
cp .github/data/fastq/test_epic_data.fastq.gz data/fastq/
```

In this case, the similarity search mode in the config file can be changed to "test", searching only one of the two strands.

## Additional Notes

The pipeline is designed to handle large sequencing datasets in parallel, so it's recommended to run it on a machine with sufficient computational resources. However, to run the pipeline on machines with less resources, it is recommended to split the fastq files or the tables in smaller chunks to prevent RAM overflow. This can be done by increase num_parts in the config file.
If any errors occur during the pipeline run, Snakemake will provide detailed logs, allowing you to debug and troubleshoot any issues. You are most welcome to create an Issue when running into problems.

## Connected 16S Analysis

We acknowledge that many users perform a regular 16S experiment additional to the epicPCR experiment. However, we have decided not to include an extra analysis  featue for this in ERMA. We would like to encourage the user to use RiboSnake (10.46471/gigabyte.132), a validated, automated, reproducible QIIME2-based pipeline implemented in Snakemake for analysing 16S rRNA gene amplicon sequencing data.

## License

This project is licensed under the MIT License.
