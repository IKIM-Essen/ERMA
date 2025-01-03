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

def filter_blast_results(input_file, output_file, min_similarity):
    df = pd.read_csv(input_file, header=0, sep=',',dtype=dtype_dict)
    
    # Filter results based on percentage identity and alignment length for ABR
    abr_data = df[
        (df['part'] == 'ABR') &
        (df['perc_identity'] > float(min_similarity) * 100) 
        #& (df['orientation'] == 'forward')
    ]
    abr_data = abr_data.groupby("query_id").apply(process_orientation_and_counts).reset_index(drop=True)

    # Filter results based on percentage identity and alignment length for 16S
    sixteen_s_data = df[
        (df['part'] == '16S')
    ]
    
    sixteen_s_data = sixteen_s_data.copy()
    sixteen_s_data["query_id"] = sixteen_s_data["query_id"].str.split(expand=True)[0]
    sixteen_s_data = sixteen_s_data.groupby("query_id").apply(process_orientation_and_counts).reset_index(drop=True)
    
    # Filter for common query IDs between abr_data and sixteen_s_data
    common_query_ids = pd.Index(abr_data['query_id']).intersection(sixteen_s_data['query_id'])
    abr_data_filtered = abr_data[abr_data['query_id'].isin(common_query_ids)]
    sixteen_s_data_filtered = sixteen_s_data[sixteen_s_data['query_id'].isin(common_query_ids)]
    
    # Concatenate the filtered ABR and 16S data
    merged_data = pd.concat([abr_data_filtered, sixteen_s_data_filtered])
    
    # Write the filtered data to output
    merged_data.to_csv(output_file, index=False)

if __name__ == "__main__":
    input_file = snakemake.input.integrated_data
    output_file = snakemake.output.filtered_data
    min_similarity = snakemake.params.min_similarity
    sys.stderr = open(snakemake.log[0], "w")
    filter_blast_results(input_file, output_file,min_similarity)
