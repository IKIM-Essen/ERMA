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

def process_card_results(card_results_path, aro_mapping_path, output_path, chunksize):
    """Process CARD results and save them to an intermediate output file."""
    blast_columns = ["query_id", "subject_id", "perc_identity", "align_length", "mismatches",
                     "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
    
    aro_df = pd.read_csv(aro_mapping_path, sep="\t")
    header_written = False

    with gzip.open(card_results_path, 'rt') as f_in, open(output_path, 'w') as f_out:
        for chunk in pd.read_csv(f_in, sep="\t", names=blast_columns, chunksize=chunksize):
            chunk["part"] = "ABR"
            chunk['ARO Accession'] = chunk['subject_id'].str.split("|", expand=True)[2]
            chunk["distance"] = chunk["q_start"] - chunk["q_end"]
            orientation_counts = chunk.groupby("query_id").apply(process_orientation_and_counts).reset_index()
            merged_chunk = chunk.merge(orientation_counts, on="query_id")
            merged_chunk = merged_chunk.merge(aro_df, on='ARO Accession', how='left')
            merged_chunk.to_csv(f_out, index=False, header=not header_written)
            header_written = True

def process_silva_results(silva_results_path, output_path, chunksize):
    """Process SILVA results and save them to an intermediate output file."""
    blast_columns = ["query_id", "subject_id", "perc_identity", "align_length", "mismatches",
                     "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
    header_written = False

    with gzip.open(silva_results_path, 'rt') as f_in, open(output_path, 'w') as f_out:
        for chunk in pd.read_csv(f_in, sep="\t", names=blast_columns, chunksize=chunksize):
            chunk["part"] = "16S"
            chunk["primaryAccession"] = chunk['subject_id'].str.split(".", expand=True)[0]
            chunk["distance"] = chunk["q_start"] - chunk["q_end"]
            orientation_counts = chunk.groupby("query_id").apply(process_orientation_and_counts).reset_index()
            merged_chunk = chunk.merge(orientation_counts, on="query_id")
            merged_chunk["genus"] = merged_chunk["subject_id"].str.split(';').str[-2]
            merged_chunk.to_csv(f_out, index=False, header=not header_written)
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
    chunksize = snakemake.params.chunksize
    sys.stderr = open(snakemake.log[0], "w")

    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_card = executor.submit(process_card_results, card_results, aro_mapping, card_output, chunksize)
        future_silva = executor.submit(process_silva_results, silva_results, silva_output, chunksize)
        
        future_card.result()
        future_silva.result()

    merge_results(card_output, silva_output, final_output)
    
    print(f"Processing complete. Final merged output saved to {final_output}")