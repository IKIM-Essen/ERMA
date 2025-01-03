# ERMA - epicPCR Resistome-Microbiome Analyzer

This pipeline is designed to process sequencing reads obtained from epicPCR experiments, linking antimicrobial resistance (AMR) genes with 16S rRNA genes. The pipeline uses **Snakemake** to manage the workflow and integrates tools for downloading necessary databases, running sequence alignments, filtering, and generating visual reports.

## Features
- Downloads and prepares SILVA and CARD databases for use in the analysis.
- Converts raw FASTQ sequencing files to FASTA format.
- Performs BLAST sequence alignments against both CARD (AMR genes) and SILVA (16S rRNA).
- Integrates BLAST results to identify linked AMR and microbial markers.
- Generates filtered results based on alignment quality and similarity thresholds.
- Produces graphical outputs such as genus distribution plots, e-value boxplots, and read positions.
- Generates an HTML report summarizing the analysis.

## Prerequisites

### Install Snakemake
The pipeline requires **Snakemake** with support for conda environments. You can install it using conda (via Miniconda or Anaconda):

```bash
conda install -c bioconda -c conda-forge snakemake
```

Alternatively, follow the official [Snakemake installation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) guide for more options.

Install Dependencies
This pipeline uses conda environments to manage its dependencies. Snakemake will automatically create and manage these environments when run with the 
`--use-conda` flag.

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
  download-path-seq: "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_LSUParc_tax_silva.fasta.gz"
  download-path-tax: "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_lsu_parc_138.2.txt.gz"

card:
  download-path: "https://card.mcmaster.ca/download/0/broadstreet-v3.3.0.tar.bz2"

chunksize: 500000 # number of lines per chunk the large tables are split
num_parts: 10 # number of chunks the fastqs are split into
max_threads: 32
```
Replace /path/to/your/project with the actual path to your local pipeline directory.
Run the Pipeline: To start the pipeline, run the following command from the base directory:

```bash
snakemake --use-conda --cores N
```
Replace N with the number of cores (threads) you want to use. Dont forget to also specify the number of threads in the config file.

## Additional Notes

The pipeline is designed to handle large sequencing datasets in parallel, so it's recommended to run it on a machine with sufficient computational resources. However, to run the pipeline on machines with less resources, it is recommended to split the fastq files or the tables in smaller chunks to prevent the RAM to overflow. Care: the lower the number of lines per chunk for the tables, the higher the chance some blast results for the same read will be lost.
If any errors occur during the pipeline run, Snakemake will provide detailed logs, allowing you to debug and troubleshoot any issues. You are most welcome to create an Issue when running into problems.

License

This project is licensed under the MIT License.
