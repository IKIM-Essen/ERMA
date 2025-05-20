import pandas as pd
import os
import sys

"""
This script filters and merges BLAST result data from ABR and 16S rRNA hits. The exact filtering steps are

- Filtering ABR and 16S entries based on user-defined percent identity thresholds
- Retaining only the highest-identity hit per query ID
- Keeping only hits with overlapping query ID in both databases

In case, one of the targets is empty, a dummy is created preventing workflow failure
"""

dtype_dict = {
    "query_id": "string",
    "perc_identity": "float",
    "align_length": "int",
    "evalue": "float",
    "part": "string",
    "genus": "string",
    "AMR Gene Family": "string",
}


def write_dummy_line(output_file):
    """Write a dummy line to the output if the input is a placeholder 16S hit"""
    print("Detected only a dummy 16S line — generating merged dummy output.")
    dummy_line = {
        "query_id": "dummy",
        "perc_identity": 0,
        "align_length": 0,
        "evalue": 0,
        "part": "16S",
        "genus": "Unclassified",
        "AMR Gene Family": "NA",
    }
    dummy_df = pd.DataFrame(dummy_line, index=[0])
    dummy_df.to_csv(output_file, index=False)
    return


def read_input_data(input_file):
    """Load relevant columns from input file with proper dtypes"""
    return pd.read_csv(input_file, sep=",", dtype=dtype_dict, usecols=dtype_dict.keys())


def filter_by_identity(df, part, min_similarity):
    """Filter BLAST result for either ABR or 16S part based on percent identity"""
    data_pre = df[df["part"] == part]
    filtered = data_pre[data_pre["perc_identity"] > min_similarity * 100]
    filtered_count = len(data_pre) - len(filtered)
    return filtered, filtered_count


def keep_max_identity_per_query(df):
    """For each query_id, keep only rows with the highest percent identity"""
    max_identities = df.groupby("query_id")["perc_identity"].max().reset_index()
    merged = df.merge(max_identities, on=["query_id", "perc_identity"])
    return merged


def clean_16s_query_ids(df):
    """Remove anything after the first whitespace in 16S query IDs"""
    df["query_id"] = df["query_id"].str.split().str[0]
    return df


def merge_parts_on_query_id(abr_data, s16_data):
    """Return only rows with query_ids present in both ABR and 16S data"""
    common_ids = pd.Index(abr_data["query_id"]).intersection(s16_data["query_id"])
    return (
        abr_data[abr_data["query_id"].isin(common_ids)],
        s16_data[s16_data["query_id"].isin(common_ids)],
    )


def write_summary(overview_table, sample, part, stats):
    """Write all filtering summary statistics to the overview file"""
    with open(overview_table, "a") as file:
        for stat_name, value in stats.items():
            file.write(f"{stat_name},{sample},{part},{value}\n")


def filter_blast_results(input_file, output_file, min_similarity, overview_table):
    """Main filtering logic for BLAST results across ABR and 16S data parts"""
    df = read_input_data(input_file)
    df_overview = pd.read_csv(
        overview_table, names=["state", "sample", "No", "total_count"]
    )

    # ABR filtering
    abr_filtered, abr_removed_identity = filter_by_identity(df, "ABR", min_similarity)
    abr_final = keep_max_identity_per_query(abr_filtered)
    abr_removed_max = len(abr_filtered) - len(abr_final)

    # 16S filtering
    s16_filtered, s16_removed_identity = filter_by_identity(df, "16S", min_similarity)
    s16_filtered = clean_16s_query_ids(s16_filtered)
    s16_final = keep_max_identity_per_query(s16_filtered)
    s16_removed_max = len(s16_filtered) - len(s16_final)

    # Handle dummy 16S result
    if len(s16_final) == 1 and s16_final.iloc[0]["query_id"] == "dummy.dummy":
        write_dummy_line(output_file)
        merge_output = df_overview.loc[
            df_overview["state"] == "merge output", "total_count"
        ].values[0]
        filtered = (
            abr_removed_identity
            + abr_removed_max
            + s16_removed_identity
            + s16_removed_max
        )
        remaining = merge_output - filtered
        sample, part = os.path.normpath(input_file).split(os.sep)[-3:-1]
        # Write summary in case of dummy
        stats = {
            "filtered min similarity ABR": "-" + str(abr_removed_identity),
            "filtered max identity ABR": "-" + str(abr_removed_max),
            "filtered min similarity 16S": "-" + str(s16_removed_identity),
            "filtered max identity 16S": "-" + str(s16_removed_max),
            "filtered query id mismatch": "-" + str(remaining),
            "filtration output": 1,
        }
        write_summary(overview_table, sample, part, stats)
        return

    # Match ABR and 16S by query_id
    abr_common, s16_common = merge_parts_on_query_id(abr_final, s16_final)
    removed_query_id_mismatch = (len(abr_final) + len(s16_final)) - (
        len(abr_common) + len(s16_common)
    )

    # Merge and write final output
    merged = pd.concat([abr_common, s16_common])
    merged.to_csv(output_file, index=False)

    # Extract sample and part from file path
    sample, part = os.path.normpath(input_file).split(os.sep)[-3:-1]

    # Write summary
    stats = {
        "filtered min similarity ABR": "-" + str(abr_removed_identity),
        "filtered max identity ABR": "-" + str(abr_removed_max),
        "filtered min similarity 16S": "-" + str(s16_removed_identity),
        "filtered max identity 16S": "-" + str(s16_removed_max),
        "filtered query id mismatch": "-" + str(removed_query_id_mismatch),
        "filtration output": len(merged),
    }
    write_summary(overview_table, sample, part, stats)


if __name__ == "__main__":
    input_file = snakemake.input.integrated_data
    overview_table = snakemake.input.overview_table
    output_file = snakemake.output.filtered_data
    min_similarity = snakemake.params.min_similarity
    sys.stderr = open(snakemake.log[0], "w")
    filter_blast_results(input_file, output_file, float(min_similarity), overview_table)
