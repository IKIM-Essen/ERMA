import pandas as pd
import os

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
    print("Detected only a dummy 16S line — generating merged dummy output.")
    dummy_line = {
        'query_id': 'dummy',
        'perc_identity': 0,
        'align_length': 0,
        'evalue': 0,
        'part': '16S',
        'genus': 'Unclassified',
        'AMR Gene Family': 'NA',
    }
    merged_data = pd.DataFrame([dummy_line])
    merged_data.to_csv(output_file, index=False)
    return  # Exit the function early

def filter_blast_results(input_file, output_file, min_similarity, overview_table):
    df = pd.read_csv(input_file, header=0, sep=',', dtype=dtype_dict, usecols=["query_id", "perc_identity", "part", "align_length","evalue","genus","AMR Gene Family"])
    
    # Filter results based on percentage identity for ABR
    abr_data_pre = df[df['part'] == 'ABR']
    abr_filtered = abr_data_pre[abr_data_pre["perc_identity"] > float(min_similarity) * 100]
    filtered_abr_identity = len(abr_data_pre) - len(abr_data)
    
    abr_summary = abr_filtered.groupby("query_id").agg(max_perc_identity=("perc_identity", "max")).reset_index()
    abr_merged = abr_filtered.merge(abr_summary,on="query_id",how="inner")
    abr_data = abr_merged[abr_merged["perc_identity"] == abr_merged["max_perc_identity"]].copy()

    filtered_abr_max = len(abr_filtered) - len(abr_data)

    # Filter results based on percentage identity for 16S
    sixteen_s_data_pre = df[df['part'] == '16S']
    sixteen_s_filtered = sixteen_s_data_pre[sixteen_s_data_pre['perc_identity'] > float(min_similarity) * 100].copy()
    sixteen_s_filtered["query_id"] = sixteen_s_filtered["query_id"].str.split().str[0]
    filtered_16s_identity = len(sixteen_s_data_pre) - len(sixteen_s_data)
    
    sixteen_s_summary = sixteen_s_filtered.groupby("query_id").agg(max_perc_identity=("perc_identity", "max")).reset_index()
    sixteen_s_merged = sixteen_s_filtered.merge(sixteen_s_summary,on="query_id",how="inner")
    sixteen_s_data = sixteen_s_merged[sixteen_s_merged["perc_identity"] == sixteen_s_merged["max_perc_identity"]].copy()
    
    filtered_16s_max = len(sixteen_s_filtered) - len(sixteen_s_data)
    print(f"filtered_max_identity_16S: {filtered_16s_max}")

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
        sample, part = os.path.normpath(input_file).split(os.sep)[-3:-1]

        # Count number of rows in the final DataFrame
        final_count = len(merged_data)

        # Write summary lines
        with open(overview_table, "a") as file:
            file.write(f"filtered_min_similarity_ABR,{sample},{part},{filtered_abr_identity}\n")
            file.write(f"filtered_max_identity_ABR,{sample},{part},{filtered_abr_max}")            
            file.write(f"filtered_min_similarity_16S,{sample},{part},{filtered_16s_identity}\n")
            file.write(f"filtered_max_identity_16S,{sample},{part},{filtered_16s_max}")            
            file.write(f"filtered_query_id_mismatch,{sample},{part},{removed_due_to_query_id}\n")
            file.write(f"filtration_output,{sample},{part},{final_count}\n")   

if __name__ == "__main__":
    input_file = snakemake.input.integrated_data
    overview_table = snakemake.input.overview_table
    output_file = snakemake.output.filtered_data
    min_similarity = snakemake.params.min_similarity
    sys.stderr = open(snakemake.log[0], "w")
    filter_blast_results(input_file, output_file,min_similarity,overview_table)
