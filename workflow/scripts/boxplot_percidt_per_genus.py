# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys


def generate_percentage_idt_per_genus(input_files, output_file):
    all_data = []  # List to hold DataFrames from all input files

    for input_file in input_files:
        df = pd.read_csv(
            input_file,
            sep=",",
            header=0,
            compression="gzip",
        )
        all_data.append(df)

    # Combine all partitions into a single DataFrame
    combined_data = pd.concat(all_data)

    # Calculate genus query counts
    genus_query_counts = (
        combined_data.groupby("genus")["query_id"].nunique().reset_index()
    )
    genus_query_counts.columns = ["genus", "unique_query_count"]

    # Keep only the top 20 genera
    top20_species = genus_query_counts.nlargest(20, "unique_query_count")

    # Filter combined_data to retain only the top 20 genera
    combined_data = combined_data[combined_data["genus"].isin(top20_species["genus"])]

    # Now filter genus_query_counts as well
    genus_query_counts = genus_query_counts[
        genus_query_counts["genus"].isin(top20_species["genus"])
    ]

    # Define order for the x-axis
    genus_order = top20_species.sort_values(by="unique_query_count", ascending=False)[
        "genus"
    ]

    # Plotting
    fig, ax1 = plt.subplots(figsize=(15, 8))
    sns.boxplot(
        x="genus",
        y="perc_identity_16S",
        data=combined_data,
        ax=ax1,
        order=genus_order,
        fliersize=0.0,
        color="dodgerblue",
    )
    ax1.set_xlabel("Bacterial Genus")
    ax1.set_ylabel("Percentage Identity (boxplot)", color="royalblue")
    ax1.set_title(
        "Boxplot of Percentage Identity and Read Counts for Each Bacterial Genus"
    )
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)

    # Add a second y-axis for unique query counts
    ax2 = ax1.twinx()
    sns.barplot(
        x="genus",
        y="unique_query_count",
        data=genus_query_counts,
        ax=ax2,
        alpha=0.2,
        color="purple",
        order=genus_order,
    )
    ax2.set_ylabel("Number of fusion reads (bar)", color="violet")

    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


if __name__ == "__main__":
    input_files = sorted(
        snakemake.input.filtered_data
    )  # List of all filtered fasta files
    output_file = snakemake.output[0]
    sys.stderr = open(snakemake.log[0], "w")
    generate_percentage_idt_per_genus(input_files, output_file)
