import pandas as pd
import os
from glob import glob

# Define necessary columns
necessary_columns = [
    "query_id",
    "part",
    "genus",
    "AMR Gene Family",
    "perc_identity",
]


def write_dummy_line(sample_name):
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
    # Separate ABR and 16S data
    abr_data = combined_data[combined_data["part"] == "ABR"]
    sixteen_s_data = combined_data[combined_data["part"] == "16S"]

    # Dummy handling
    if sixteen_s_data.empty or abr_data.empty:
        return write_dummy_line(sample_name)

    # Drop duplicates to ensure 1:1 merge
    abr_unique = abr_data[["query_id", "AMR Gene Family"]].drop_duplicates()
    sixteen_unique = sixteen_s_data[["query_id", "genus"]].drop_duplicates()

    merged = pd.merge(abr_unique, sixteen_unique, on="query_id", how="inner")
    merged["sample"] = sample_name

    # Count per AMR Gene Family and genus
    genus_counts = (
        merged.groupby(["sample", "AMR Gene Family", "genus"])
        .size()
        .reset_index(name="genus_count")
    )

    # Total genus count per AMR Gene Family
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
    data_frames = []
    for file in file_list:
        try:
            df = pd.read_csv(file, usecols=necessary_columns, compression="gzip")
            data_frames.append(df)
        except Exception as e:
            print(f"Skipping {file}: {e}")
    if data_frames:
        merged_df = pd.concat(data_frames, ignore_index=True)
    else:
        merged_df = pd.DataFrame(columns=necessary_columns)
    return merged_df


def export_genera_abundance(input_files, output_html):
    # Group input files by sample
    sample_to_files = {}
    for file in input_files:
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
    final_df.to_csv(output_csv, index=False)


if __name__ == "__main__":
    input_files = list(snakemake.input.filtered_data)
    output_csv = snakemake.output.csv
    sys.stderr = open(snakemake.log[0], "w")
    export_genera_abundance(input_files, output_csv)
