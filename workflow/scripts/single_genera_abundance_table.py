# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


import pandas as pd
import os, sys

"""
This script processes epicPCR data (ABR + 16S) distingly for one sample and its part
to compute genus-level total and relative abundance per AMR Gene Family in a table.
This table is later be imported as sample-wise information in the report.
"""


def write_dummy_line(sample_name):
    # Create a dummy row when either ABR or 16S data is missing for a sample
    # Returns Single-row dataframe with placeholder values
    dummy_line = {
        "sample": sample_name,
        "AMR Gene Family": "NA",
        "genus": "NA",
        "genus_count": 0,
        "total_count": 0,
        "relative_genus_count": 0,
    }
    merged_data = pd.DataFrame([dummy_line])
    return merged_data


def process_combined_data(combined_data, sample_name):
    # Separate ABR and 16S data for merging by query_id
    df = combined_data
    if df.iloc[0]["query_id"] == "dummy" or df.empty:
        return write_dummy_line(sample_name)

    # Add the sample name
    df["sample"] = sample_name

    # Calculate genus counts per AMR Gene Family and genus for the sample
    genus_counts = (
        df.groupby(["sample", "AMR Gene Family", "genus"])
        .size()
        .reset_index(name="genus_count")
    )

    # Calculate total genus count per AMR Gene Family within each sample
    total_counts = (
        genus_counts.groupby(["sample", "AMR Gene Family"])["genus_count"]
        .sum()
        .reset_index(name="total_count")
    )

    # Merge to get total counts for each genus entry and calculate relative counts
    genus_counts = pd.merge(
        genus_counts, total_counts, on=["sample", "AMR Gene Family"], how="left"
    )
    genus_counts["relative_genus_count"] = round(
        genus_counts["genus_count"] / genus_counts["total_count"], 4
    )

    return genus_counts


def export_genera_abundance(input_files, sample_name, parts, output_path):
    sample_input_files = [f for f in input_files if f"/{sample_name}/" in f]
    # Load and combine all parts for the current sample
    part_dfs = []
    for part in parts:
        matching_files = [f for f in sample_input_files if f"/{part}/" in f]
        if not matching_files:
            continue
        input_file = matching_files[0]
        df = pd.read_csv(input_file, sep=",", header=0, compression="gzip")
        part_dfs.append(df)

    if not part_dfs:
        print(f"No valid parts found for sample: {sample_name}")
        dummy_df = write_dummy_line(sample_name)
        dummy_df.to_csv(output_path, index=False)
        return

    full_sample_df = pd.concat(part_dfs, ignore_index=True)
    processed_data = process_combined_data(full_sample_df, sample_name)

    processed_data = processed_data.sort_values(
        by=["total_count", "genus_count"], ascending=False
    )

    # Write to HTML
    processed_data.to_html(output_path, index=False)


if __name__ == "__main__":
    input_file = snakemake.input.filtered_data
    output_path = snakemake.output[0]
    sample_name = snakemake.params.sample_name
    parts = snakemake.params.parts
    sys.stderr = open(snakemake.log[0], "w")
    export_genera_abundance(input_file, sample_name, parts, output_path)
