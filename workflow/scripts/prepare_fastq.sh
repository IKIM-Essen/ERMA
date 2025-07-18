#!/bin/bash
# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.

: '
This script is meant to merge paired-end Illumina reads to one file.
The script loads every R1/R2_fastq.gz file and merges them together.
It is meant to use before the ERMA pipeline as preparation.
'

# Define paths
FASTQ_DIR="data/fastq"
MERGED_DIR="${FASTQ_DIR}/input"
SM_DIR=$(pwd)

# Create the merged directory if it doesn't exist
mkdir -p "$MERGED_DIR"

# Copy all fastq.gz files into the subfolder
mv ${FASTQ_DIR}/*.fastq.gz "$MERGED_DIR"

# Move into the merged directory
cd "$MERGED_DIR" || exit

# Merge all R1/R2 pairs with NGmerge
for r1 in *_R1_*.fastq.gz; do
    # Define the corresponding R2 file
    r2="${r1/_R1_/_R2_}"
    
    # Extract the sample name (remove _R1_001)
    sample_name=$(echo "$r1" | sed 's/_R1_.*//' | sed 's/_L001//')
    echo "Preparing: $sample_name"
    # Run NGmerge
    echo "Merging $r1 and $r2 into ${sample_name}.fastq.gz..."
    NGmerge -1 "$r1" -2 "$r2" -o "../${sample_name}.fastq.gz"

    # Optional: Remove original files after merging
    # rm "$r1" "$r2"
done

# Move back to project root
cd "$SM_DIR" || exit
touch "$FASTQ_DIR/merged_done.txt"
