import pandas as pd
import sys

dtype_dict = {
    "query_id": "string",
    "subject_id": "string",
    "perc_identity": "float",
    "align_length": "int",
    "mismatches": "int",
    "gap_opens": "int",
    "q_start": "int",
    "q_end": "int",
    "s_start": "int",
    "s_end": "int",
    "evalue": "float",
    "bit_score": "float",
    "part": "string",
    "primaryAccession": "string",
    "acc_start": "float",
    "acc_stop": "float",
    "distance": "int",
    "orientation": "string",
    "start": "float",
    "stop": "float",
    "organism_name": "string",
    "taxid": "float",
    "ARO Accession": "string",
    "CVTERM ID": "float",
    "Model Sequence ID": "float",
    "Model ID": "float",
    "Model Name": "string",
    "ARO Name": "string",
    "Protein Accession": "string",
    "DNA Accession": "string",
    "AMR Gene Family": "string",
    "Drug Class": "string",
    "Resistance Mechanism": "string",
    "CARD Short Name": "string"
}

def write_dummy_line(output_file):
    print("Detected only a dummy 16S line — generating merged dummy output.")
    dummy_line = {
        'query_id': 'dummy',
        'subject_id': 'NA',
        'perc_identity': 0,
        'align_length': 0,
        'mismatches': 0,
        'gap_opens': 0,
        'q_start': 0,
        'q_end': 0,
        's_start': 0,
        's_end': 0,
        'evalue': 0,
        'bit_score': 0,
        'part': '16S',
        'primaryAccession': 'NA',
        'distance': 0,
        'orientation': 'NA',
        'genus': 'Unclassified',
        'ARO Accession': 'NA',
        'CVTERM ID': 'NA',
        'Model Sequence ID': 'NA',
        'Model ID': 'NA',
        'Model Name': 'NA',
        'ARO Name': 'NA',
        'Protein Accession': 'NA',
        'DNA Accession': 'NA',
        'AMR Gene Family': 'NA',
        'Drug Class': 'NA',
        'Resistance Mechanism': 'NA',
        'CARD Short Name': 'NA',
        'most_common_q_start': 0,
        'most_common_q_end': 0
    }
    merged_data = pd.DataFrame([dummy_line])
    merged_data.to_csv(output_file, index=False)
    return  # Exit the function early

def filter_group(group):
    most_common_q_start = group["most_common_q_start"].dropna().iloc[0]
    abr_part = group[group["part"] == "ABR"]
    filtered_16s_part = group[
        (group["part"] == "16S") & 
        (group["q_start"].between(most_common_q_start - 30, most_common_q_start + 30))
    ]
    if filtered_16s_part.empty:
        return pd.DataFrame()     
    return pd.concat([abr_part, filtered_16s_part])

def process_orientation_and_counts(group):
    # Get the most common start and end positions
    most_common_q_start = group["q_start"].mode().iloc[0]
    most_common_q_end = group["q_end"].mode().iloc[0]
    
    # Filter for the maximum perc_identity rows within this group
    max_perc_identity = group["perc_identity"].max()
    filtered_group = group[group["perc_identity"] == max_perc_identity]

    # Add the start/end positions as columns for each row in the filtered group
    filtered_group = filtered_group.assign(
        most_common_q_start=most_common_q_start,
        most_common_q_end=most_common_q_end
    )
    return filtered_group

def filter_blast_results(input_file, output_file, min_similarity, overview_table):
    df = pd.read_csv(input_file, header=0, sep=',', dtype=dtype_dict)
    
    total_rows = len(df)
    
    # Filter results based on percentage identity for ABR
    abr_data_pre = df[df['part'] == 'ABR']
    abr_data = abr_data_pre[abr_data_pre['perc_identity'] > float(min_similarity) * 100]
    filtered_abr_identity = len(abr_data_pre) - len(abr_data)
    
    abr_data = abr_data.groupby("query_id").apply(process_orientation_and_counts).reset_index(drop=True)

    # Filter results based on percentage identity for 16S
    sixteen_s_data_pre = df[df['part'] == '16S']
    sixteen_s_data = sixteen_s_data_pre[sixteen_s_data_pre['perc_identity'] > float(min_similarity) * 100]
    filtered_16s_identity = len(sixteen_s_data_pre) - len(sixteen_s_data)
    
    sixteen_s_data = sixteen_s_data.copy()
    sixteen_s_data["query_id"] = sixteen_s_data["query_id"].str.split(expand=True)[0]
    sixteen_s_data = sixteen_s_data.groupby("query_id").apply(process_orientation_and_counts).reset_index(drop=True)
    
    if len(sixteen_s_data) == 1 and sixteen_s_data.iloc[0]["query_id"] == "dummy.dummy":
        write_dummy_line(output_file)
    else:
        # Filter for common query IDs
        common_query_ids = pd.Index(abr_data['query_id']).intersection(sixteen_s_data['query_id'])
        abr_data_filtered = abr_data[abr_data['query_id'].isin(common_query_ids)]
        sixteen_s_data_filtered = sixteen_s_data[sixteen_s_data['query_id'].isin(common_query_ids)]

        removed_due_to_query_id = (len(abr_data) + len(sixteen_s_data)) - (len(abr_data_filtered) + len(sixteen_s_data_filtered))

        # Concatenate the filtered ABR and 16S data
        merged_data = pd.concat([abr_data_filtered, sixteen_s_data_filtered])
        merged_data.to_csv(output_file, index=False)

        # Extract sample and part from file path
        path_parts = input_file.split("/")
        sample = path_parts[1]
        part = path_parts[2]

        # Count number of rows in the final DataFrame
        final_count = len(merged_data)

        # Write summary lines
        with open(overview_table, "a") as file:
            file.write(f"filtered_min_similarity_ABR,{sample},{part},{filtered_abr_identity}\n")
            file.write(f"filtered_min_similarity_16S,{sample},{part},{filtered_16s_identity}\n")
            file.write(f"filtered_query_id_mismatch,{sample},{part},{removed_due_to_query_id}\n")
            file.write(f"filtration_output,{sample},{part},{final_count}\n")   

if __name__ == "__main__":
    input_file = snakemake.input.integrated_data
    overview_table = snakemake.input.overview_table
    output_file = snakemake.output.filtered_data
    min_similarity = snakemake.params.min_similarity
    sys.stderr = open(snakemake.log[0], "w")
    filter_blast_results(input_file, output_file,min_similarity,overview_table)
