# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

"""
This script takes a list of all filtered fasta files, combines e-value information 
across samples, and visualizes the distribution of e-values using boxplots split 
by part (ABR/16S) and sample.
"""


def read_and_process_partitioned_data(partition_files, sample):
    """
    Read and process filtered data files for a given sample.

    Args:
        partition_files (list of str): Paths to CSV files corresponding to sample parts.
        sample (str): Sample name.

    Returns:
        pd.DataFrame or None: Combined DataFrame of filtered entries or None if files missing.
    """
    data_frames = []
    sample_name = sample

    for part_file in partition_files:
        if os.path.exists(part_file):
            fields = ["query_id", "evalue", "part"]
            df = pd.read_csv(
                part_file, header=0, sep=",", usecols=fields, compression="gzip"
            )
            df = df.drop_duplicates()
            df["sample"] = sample_name
            data_frames.append(df)
        else:
            print(f"File {part_file} does not exist.")

    if data_frames:
        return pd.concat(data_frames)
    else:
        return None


def plot_boxplots(data, output_file):
    """
    Generate and save boxplots of e-values across samples and parts (ABR vs. 16S).

    Args:
        data (pd.DataFrame): Combined dataframe containing 'sample', 'evalue', and 'part'.
        output_file (str): Path to save the resulting plot.
    """
    plt.figure(figsize=(15, 10))
    flierprops = dict(markerfacecolor="0.75", markersize=2, linestyle="none")
    sns.boxplot(x="sample", y="evalue", hue="part", data=data, flierprops=flierprops)
    plt.yscale("log")
    plt.title("Boxplot of e-values for ABR and 16S parts across samples -Filtered-")
    plt.xlabel("Sample")
    plt.ylabel("E-value (log scale)")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


def main(filtered_fasta_files, sample_names, output_file):
    """Main function to process partitioned files for each sample and generate the plot."""
    all_data = []

    # Loop over each sample's partitioned CSV files
    for sample in sample_names:
        data = read_and_process_partitioned_data(
            [file for file in filtered_fasta_files if str(sample) in file], sample
        )
        if data is not None:
            all_data.append(data)

    if all_data:
        combined_data = pd.concat(all_data)
        plot_boxplots(combined_data, output_file)
    else:
        print("No data found.")


if __name__ == "__main__":
    filtered_fasta_files = sorted(
        snakemake.input.filtered_data
    )  # List of all filtered fasta files files
    output_file = snakemake.output[0]  # Path to save the output plot
    sample_name = sorted(snakemake.params.sample_name)  # Minimum similarity filter
    sys.stderr = open(snakemake.log[0], "w")
    main(filtered_fasta_files, sample_name, output_file)
