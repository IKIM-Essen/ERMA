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
This script takes a list of all filtered fasta files, combines param information 
across samples, and visualizes the distribution of param using boxplots split 
by part (ABR/16S) and sample.
"""

PRETTY_LABELS = {
    "align_length": "Alignment length",
    "perc_identity": "Percentage identity",
    "evalue": "E-value",
}


def read_and_process_partitioned_data(partition_files, sample, param):
    """Read and process partitioned files for a single sample."""
    data_frames = []
    sample_name = sample
    param = param
    for part_file in partition_files:
        if os.path.exists(part_file):
            df = pd.read_csv(part_file, header=0, sep=",")
            long_df = pd.melt(
                df,
                id_vars=["query_id"],
                value_vars=[param + "_ABR", param + "_16S"],
                var_name="part",
                value_name=param,
            )

            # Normalize part labels
            long_df["part"] = long_df["part"].str.replace(param + "_", "")
            long_df["sample"] = sample_name
            data_frames.append(long_df)

    if data_frames:
        return pd.concat(data_frames)
    else:
        return None


def plot_boxplots(data, param, output_file):
    """
    Generate and save boxplots of param across samples and parts (ABR vs. 16S).

    Args:
        data (pd.DataFrame): Combined dataframe containing 'sample', param, and 'part'.
        output_file (str): Path to save the resulting plot.
    """
    plt.figure(figsize=(15, 10))
    flierprops = dict(markerfacecolor="0.75", markersize=2, linestyle="none")
    sns.boxplot(x="sample", y=param, hue="part", data=data, flierprops=flierprops)
    plt.yscale("log")
    plt.title(
        f"Boxplot of {PRETTY_LABELS[param]} for ABR and 16S parts across samples -Filtered-"
    )
    plt.xlabel("Sample")
    plt.ylabel(f"{PRETTY_LABELS[param]}")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


def main(filtered_fasta_files, sample_names, param, output_file):
    """Main function to process partitioned files for each sample and generate the plot."""
    all_data = []

    # Loop over each sample's partitioned CSV files
    for sample in sample_names:
        data = read_and_process_partitioned_data(
            [file for file in filtered_fasta_files if str(sample) in file],
            sample,
            param,
        )
        data = data[data[param] > 0]
        if data is not None:
            all_data.append(data)

    if all_data:
        combined_data = pd.concat(all_data)
        plot_boxplots(combined_data, param, output_file)
    else:
        print("No data found.")


if __name__ == "__main__":
    filtered_fasta_files = sorted(
        snakemake.input.filtered_data
    )  # List of all filtered fasta files files
    output_file = snakemake.output[0]  # Path to save the output boxplot
    sample_name = sorted(snakemake.params.sample_name)  # Minimum similarity filter
    sys.stderr = open(snakemake.log[0], "w")
    param = "evalue"
    main(filtered_fasta_files, sample_name, param, output_file)
