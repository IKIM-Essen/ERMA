import pandas as pd
import gzip
import concurrent.futures
import sys
import os

def process_orientation_and_counts(group):
    orientation = "mixed"
    if (group["distance"] < 0).all():
        orientation = "reverse"
    if (group["distance"] >= 0).all():
        orientation = "forward"
    return pd.Series({"orientation": orientation})

def process_card_results(card_results_path, aro_mapping_path, output_path):
    """Process CARD results and save them to an intermediate output file."""
    blast_columns = ["query_id", "subject_id", "perc_identity", "align_length", "mismatches",
                     "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
    
    aro_df = pd.read_csv(aro_mapping_path, sep="\t")

    with open(card_results_path, 'rt') as f_in, open(output_path, 'w') as f_out:
        card_df = pd.read_csv(f_in, sep="\t", names=blast_columns)
        card_df["part"] = "ABR"
        card_df['ARO Accession'] = card_df['subject_id'].str.split("|", expand=True)[2]
        card_df["distance"] = card_df["q_start"] - card_df["q_end"]
        orientation_counts = card_df.groupby("query_id").apply(process_orientation_and_counts).reset_index()
        merged_df = card_df.merge(orientation_counts, on="query_id")
        merged_df = merged_df.merge(aro_df, on='ARO Accession', how='left')
        merged_df.to_csv(f_out, index=False)

def process_silva_results(silva_results_path, output_path):
    """Process SILVA results and save them to an intermediate output file."""
    blast_columns = ["query_id", "subject_id", "perc_identity", "align_length", "mismatches",
                     "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
    header_written = False

    with open(silva_results_path, 'rt') as f_in, open(output_path, 'w') as f_out:
        silva_df = pd.read_csv(f_in, sep="\t", names=blast_columns)
        silva_df["part"] = "16S"
        silva_df["primaryAccession"] = silva_df['subject_id'].str.split(".", expand=True)[0]
        silva_df["distance"] = silva_df["q_start"] - silva_df["q_end"]
        orientation_counts = silva_df.groupby("query_id").apply(process_orientation_and_counts).reset_index()
        merged_df = silva_df.merge(orientation_counts, on="query_id")
        merged_df["genus"] = merged_df["subject_id"].str.split(';').str[-2]
        merged_df.to_csv(f_out, index=False, header=not header_written)
        header_written = True  # Ensure header is only written once


def merge_results(card_output, silva_output, final_output):
    """Merge processed CARD and SILVA results into one final output file."""
    card_df = pd.read_csv(card_output)
    silva_df = pd.read_csv(silva_output)
    
    combined_df = pd.concat([silva_df,card_df])
    
    combined_df.to_csv(final_output, index=False)

if __name__ == "__main__":
    card_results = snakemake.input.card_results
    silva_results = snakemake.input.silva_results
    aro_mapping = snakemake.input.aro_mapping
    card_output = snakemake.output.intermed_card_results
    silva_output = snakemake.output.intermed_silva_results
    final_output = snakemake.output.integrated_data
    sys.stderr = open(snakemake.log[0], "w")

    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_card = executor.submit(process_card_results, card_results, aro_mapping, card_output)
        future_silva = executor.submit(process_silva_results, silva_results, silva_output)
        
        future_card.result()
        future_silva.result()

    merge_results(card_output, silva_output, final_output)
    
    print(f"Processing complete. Final merged output saved to {final_output}")