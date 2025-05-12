import pandas as pd
import concurrent.futures
import os, sys

"""
This script processes BLAST results from CARD (AMR gene hits) and SILVA 
(16S rRNA gene hits), infers read orientation based on alignment coordinates, 
merges metadata from the CARD ARO mapping file, and outputs cleaned, 
standardized intermediate files. It then combines both parts into a 
final merged output and logs the result in an overview table. 
The script uses concurrent execution to process CARD and SILVA results 
in parallel and ensures robustness by inserting dummy data for empty input files.
"""


def write_dummy_line(output_file, part):
    """Write a dummy line to ensure compatibility with downstream analysis"""
    if part == "ABR":
        dummy_data = pd.read_csv(snakemake.input.dummy_ABR)
    elif part == "16S":
        dummy_data = pd.read_csv(snakemake.input.dummy_16S)
    pd.DataFrame(dummy_data).to_csv(output_file, index=False)


def process_card_results(
    card_results_path, aro_mapping_path, blast_columns, output_path
):
    """Process CARD results and save them to an intermediate output file"""

    aro_df = pd.read_csv(aro_mapping_path, sep="\t")

    try:
        with open(card_results_path, "rt") as f_in, open(output_path, "w") as f_out:
            if os.stat(card_results_path).st_size == 0:
                print(f"Warning: {card_results_path} is empty. Writing dummy line.")
                write_dummy_line(output_path, "ABR")
                return

            card_df = pd.read_csv(f_in, sep="\t", names=blast_columns)
            card_df["part"] = "ABR"
            # Extract ARO accession (formatted like: ARO|...|ACCESSION|...)
            card_df["ARO Accession"] = card_df["subject_id"].str.split(
                "|", expand=True
            )[2]
            card_df["distance"] = card_df["q_start"] - card_df["q_end"]
            orientation_counts = (
                card_df.groupby("query_id")["distance"]
                .agg(
                    lambda x: (
                        "reverse"
                        if (x < 0).all()
                        else "forward" if (x >= 0).all() else "mixed"
                    )
                )
                .rename("orientation")
                .reset_index()
            )
            merged_df = card_df.merge(orientation_counts, on="query_id")
            merged_df = merged_df.merge(aro_df, on="ARO Accession", how="left")
            merged_df.to_csv(f_out, index=False)
    except Exception as e:
        print(f"Error processing {card_results_path}: {e}")
        write_dummy_line(output_path, "ABR")


def process_silva_results(silva_results_path, blast_columns, output_path):
    """Process SILVA results and save them to an intermediate output file."""

    try:
        with open(silva_results_path, "rt") as f_in, open(output_path, "w") as f_out:
            if os.stat(silva_results_path).st_size == 0:
                print(f"Warning: {silva_results_path} is empty. Writing dummy line.")
                write_dummy_line(output_path, "16S")
                return

            silva_df = pd.read_csv(f_in, sep="\t", names=blast_columns)
            silva_df["part"] = "16S"
            # Extract the primary accession (before '.') from SILVA subject_id
            silva_df["primaryAccession"] = silva_df["subject_id"].str.split(
                ".", expand=True
            )[0]
            silva_df["distance"] = silva_df["q_start"] - silva_df["q_end"]
            orientation_counts = (
                silva_df.groupby("query_id")["distance"]
                .agg(
                    lambda x: (
                        "reverse"
                        if (x < 0).all()
                        else "forward" if (x >= 0).all() else "mixed"
                    )
                )
                .rename("orientation")
                .reset_index()
            )
            merged_df = silva_df.merge(orientation_counts, on="query_id")
            merged_df["genus"] = merged_df["subject_id"].str.split(";").str[-2]
            merged_df.to_csv(f_out, index=False)
    except Exception as e:
        print(f"Error processing {silva_results_path}: {e}")
        write_dummy_line(output_path, "16S")


def merge_results(card_output, silva_output, final_output, overview_table):
    """Merge processed CARD and SILVA results into one final output file and update overview"""
    card_df = pd.read_csv(card_output)
    silva_df = pd.read_csv(silva_output)

    combined_df = pd.concat([silva_df, card_df])
    combined_df.to_csv(final_output, index=False)

    # Extract sample and part from card_output path
    path_parts = card_output.split("/")
    sample = path_parts[1]
    part = path_parts[2]

    # Count number of rows in the combined DataFrame
    count = len(combined_df)

    with open(overview_table, "a") as file:
        line = f"integration_output,{sample},{part},{count}\n"
        file.write(line)


if __name__ == "__main__":
    card_results = snakemake.input.card_results
    silva_results = snakemake.input.silva_results
    aro_mapping = snakemake.input.aro_mapping
    overview_table = snakemake.input.overview_table
    card_output = snakemake.output.intermed_card_results
    silva_output = snakemake.output.intermed_silva_results
    final_output = snakemake.output.integrated_data
    sys.stderr = open(snakemake.log[0], "w")

    blast_columns = [
        "query_id",
        "subject_id",
        "perc_identity",
        "align_length",
        "mismatches",
        "gap_opens",
        "q_start",
        "q_end",
        "s_start",
        "s_end",
        "evalue",
        "bit_score",
    ]

    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_card = executor.submit(
            process_card_results, card_results, aro_mapping, blast_columns, card_output
        )
        future_silva = executor.submit(
            process_silva_results, silva_results, blast_columns, silva_output
        )

        future_card.result()
        future_silva.result()

    merge_results(card_output, silva_output, final_output, overview_table)
