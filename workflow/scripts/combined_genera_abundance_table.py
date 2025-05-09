import pandas as pd
import os, sys

"""
This script processes epicPCR data (ABR + 16S) combined over all samples and all parts
to compute genus-level total and relative abundance per AMR Gene Family in a table.
This table is later used to create the main result bubble plot.
"""

# Necessary columns to load in each dataframe
necessary_columns = [
    "query_id",
    "part",
    "genus",
    "AMR Gene Family",
    "perc_identity",
]


def write_dummy_line(sample_name):
    """ Create a dummy row for missing input and returns placeholder """
    dummy_line = {
        "sample": sample_name,
        "AMR Gene Family": "NA",
        "genus": "NA",
        "genus_count": 0,
        "total_genus_count": 0,
        "relative_genus_count": 0,
    }
    return pd.DataFrame([dummy_line])


def process_combined_data(combined_data, sample_name):
    """ Separate ABR and 16S data for merging by query_id """
    abr_data = combined_data[combined_data["part"] == "ABR"]
    sixteen_s_data = combined_data[combined_data["part"] == "16S"]

    # Dummy handling
    if sixteen_s_data.empty or abr_data.empty:
        return write_dummy_line(sample_name)

    # Prepare to merge only unique hits
    abr_unique = abr_data[["query_id", "AMR Gene Family"]].drop_duplicates()
    sixteen_unique = sixteen_s_data[["query_id", "genus"]].drop_duplicates()

    merged = pd.merge(abr_unique, sixteen_unique, on="query_id", how="inner")
    merged["sample"] = sample_name

    # Calculate genus counts per AMR Gene Family and genus for the sample
    genus_counts = (
        merged.groupby(["sample", "AMR Gene Family", "genus"])
        .size()
        .reset_index(name="genus_count")
    )

    # Calculate total genus count per AMR Gene Family within each sample
    total_counts = (
        genus_counts.groupby(["sample", "AMR Gene Family"])["genus_count"]
        .sum()
        .reset_index(name="total_genus_count")
    )

    # Join and calculate relative abundance
    result = pd.merge(genus_counts, total_counts, on=["sample", "AMR Gene Family"])
    result["relative_genus_count"] = round(
        result["genus_count"] / result["total_genus_count"], 4
    )

    return result


def load_and_merge_parts(file_list):
    """ Load and merges dataframes over all samples """
    data_frames = []
    for file in file_list:
        try:
            df = pd.read_csv(file, usecols=necessary_columns, compression="gzip")
            data_frames.append(df)
        except Exception as e:
            print(f"Skipping file due to read error [{file}]: {repr(e)}")
    if data_frames:
        merged_df = pd.concat(data_frames, ignore_index=True)
    else:
        merged_df = pd.DataFrame(columns=necessary_columns)
    return merged_df


def export_genera_abundance(input_files, output_path):
    """ Group input files by sample """
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
    final_df = final_df.sort_values(by=["genus_count"], ascending=False)

    # Export the final aggregated data to a CSV file
    final_df.to_csv(output_path, index=False)


if __name__ == "__main__":
    input_files = list(snakemake.input.filtered_data)
    output_path = snakemake.output.csv
    sys.stderr = open(snakemake.log[0], "w")
    export_genera_abundance(input_files, output_path)
