# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


import pandas as pd
import os, sys

"""
This script processes epicPCR data (ABR + 16S) combined over all samples and all parts
to compute genus-level total and relative abundance per AMR Gene Family in a table.
This table is later used to create the main result bubble plot.
"""

# Necessary columns to load in each dataframe


def write_dummy_line(sample_name):
    """Create a dummy row for missing input and returns placeholder"""
    dummy_line = {
        "sample": sample_name,
        "AMR Gene Family": "NA",
        "genus": "NA",
        "genus_count": 0,
        "total_count": 0,
        "relative_genus_count": 0,
    }
    return pd.DataFrame([dummy_line])


def process_combined_data(combined_data, sample_name):
    combined_data["sample"] = sample_name

    # Count genus occurrences per AMR Gene Family
    genus_counts = (
        combined_data.groupby(["sample", "AMR Gene Family", "genus"])
        .size()
        .reset_index(name="genus_count")
    )

    # Calculate total genus count per AMR Gene Family within each sample
    total_counts = (
        genus_counts.groupby(["sample", "AMR Gene Family"])["genus_count"]
        .sum()
        .reset_index(name="total_count")
    )

    # Join and calculate relative abundance
    result = pd.merge(genus_counts, total_counts, on=["sample", "AMR Gene Family"])
    result["relative_genus_count"] = round(
        result["genus_count"] / result["total_count"], 4
    )

    return result


def load_and_merge_parts(file_list):
    """Load and merges dataframes over all samples"""
    data_frames = []
    for file in file_list:
        try:
            df = pd.read_csv(file, compression="gzip")
            data_frames.append(df)
        except Exception as e:
            print(f"Skipping file due to read error [{file}]: {repr(e)}")
    if data_frames:
        merged_df = pd.concat(data_frames, ignore_index=True)
    else:
        merged_df = pd.DataFrame()
    return merged_df


def export_genera_abundance(input_files, output_path):
    """Group input files by sample"""
    sample_to_files = {}
    for file in input_files:
        # Extract sample name from the file path, assuming 3rd-to-last split is the sample name
        sample = os.path.normpath(file).split(os.sep)[-3]
        sample_to_files.setdefault(sample, []).append(file)

    all_data = []

    for sample_name, files in sample_to_files.items():
        merged_data = load_and_merge_parts(files)
        sample_data = process_combined_data(merged_data, sample_name)
        all_data.append(sample_data)

    final_df = pd.concat(all_data, ignore_index=True)
    final_df = final_df.sort_values(
        by=["sample", "AMR Gene Family", "genus_count"], ascending=False
    )

    # Export the final aggregated data to a CSV file
    final_df.to_csv(output_path, index=False)


if __name__ == "__main__":
    input_files = list(snakemake.input.filtered_data)
    output_path = snakemake.output.csv
    sys.stderr = open(snakemake.log[0], "w")
    export_genera_abundance(input_files, output_path)
