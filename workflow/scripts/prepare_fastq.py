# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


import os
import sys
import pandas as pd
import subprocess
from pathlib import Path

"""
This script is meant to demultiplex a finished ONT run.
For this, a csv is used as input containing the barcode - sample relation.
The script loads every fastq.gz file for every barcode, merges and renames them.
It is meant to use before the ERMA pipeline as preparation.
"""

# Inputs from Snakemake
source_dir = snakemake.params.run_path
barcode_csv_path = snakemake.params.sample_name_path
target_length = int(snakemake.params.target_fragment_length)
filter_intervall = float(snakemake.params.filter_intervall)
output_dir = snakemake.params.output_dir

# Set filter range around the target length
min_len = int(target_length * (1 - filter_intervall))
max_len = int(target_length * (1 + filter_intervall))

Path(output_dir).mkdir(parents=True, exist_ok=True)

# Load barcode -> sample_name mapping
barcode_df = pd.read_csv(barcode_csv_path, dtype={"barcode": str})
barcode_dict = dict(zip(barcode_df["barcode"], barcode_df["sample_name"]))

for barcode, sample_name in barcode_dict.items():
    bc_dir = os.path.join(source_dir, f"barcode{barcode.zfill(2)}")
    if not os.path.exists(bc_dir):
        print(f"Skipping missing: {bc_dir}")
        continue

    # Get all .fastq.gz files
    fastq_files = list(Path(bc_dir).rglob("*.fastq.gz"))
    if not fastq_files:
        print(f"No fastq.gz files in {bc_dir}")
        continue

    output_path = os.path.join(output_dir, f"{sample_name}.fastq.gz")

    # Build seqkit command
    cat_command = "zcat " + " ".join(str(f) for f in fastq_files)
    seqkit_command = f"seqkit seq --min-len {min_len} --max-len {max_len}"

    # Run the pipeline
    with open(output_path, "wb") as out_f:
        cat_proc = subprocess.Popen(cat_command, shell=True, stdout=subprocess.PIPE)
        seqkit_proc = subprocess.Popen(
            seqkit_command, shell=True, stdin=cat_proc.stdout, stdout=subprocess.PIPE
        )
        gzip_proc = subprocess.Popen(["gzip"], stdin=seqkit_proc.stdout, stdout=out_f)
        gzip_proc.communicate()

    print(f"{sample_name}: written filtered reads to {output_path}")

# Touch sentinel file
with open(snakemake.output.sentinel, "w") as f:
    f.write("done")
