import pandas as pd
import altair as alt

# Necessary columns to load from each CSV file
necessary_columns = [
    "query_id",
    "part",
    "genus",
    "AMR Gene Family",
    "perc_identity",
]

def process_combined_data(combined_data, sample_name):
    # Separate ABR and 16S data for merging by query_id
    abr_data = combined_data[combined_data["part"] == "ABR"]
    sixteen_s_data = combined_data[combined_data["part"] == "16S"]
    
    # Prepare to merge only unique hits
    unique_abr_data = abr_data[['query_id', 'AMR Gene Family']].drop_duplicates()
    unique_sixteen_s_data = sixteen_s_data[['query_id', 'genus']].drop_duplicates()
    
    # Merge on query_id to associate AMR Gene Family with genus information from 16S data
    merged_data = pd.merge(
        unique_abr_data[['query_id', 'AMR Gene Family']], 
        unique_sixteen_s_data[['query_id', 'genus']],
        on='query_id', 
        how='inner'
    )
    
    # Add the sample name
    merged_data['sample'] = sample_name
    
    # Calculate genus counts per AMR Gene Family and genus for the sample
    genus_counts = merged_data.groupby(['sample', 'AMR Gene Family', 'genus']).size().reset_index(name='genus_count')
    
    # Calculate total genus count per AMR Gene Family within each sample
    total_counts = genus_counts.groupby(['sample', 'AMR Gene Family'])['genus_count'].sum().reset_index(name='total_genus_count')
    
    # Merge to get total counts for each genus entry and calculate relative counts
    genus_counts = pd.merge(genus_counts, total_counts, on=['sample', 'AMR Gene Family'], how='left')
    genus_counts['relative_genus_count'] = round(genus_counts['genus_count'] / genus_counts['total_genus_count'],4)
    
    return genus_counts

def combine_blast_data(input_files, sample_name):
    # Load and combine data from all parts for the given sample
    all_data = [pd.read_csv(input_file, sep=",", usecols=necessary_columns, header=0, compression='gzip') for input_file in input_files]
    combined_data = pd.concat(all_data, ignore_index=True)
    
    # Process combined data to get genus counts and relative values
    genus_counts = process_combined_data(combined_data, sample_name)
    return genus_counts

def export_genera_abundance(input_files, sample_names, output_csv):
    all_samples_data = pd.DataFrame()
    
    # Process each sample’s files to build the final DataFrame
    for sample_name in sample_names:
        sample_files = [f for f in input_files if f"/{sample_name}/" in f]
        if sample_files:
            sample_data = combine_blast_data(sample_files, sample_name)
            all_samples_data = pd.concat([all_samples_data, sample_data], ignore_index=True)
    
    # Export the final aggregated data to a CSV file
    all_samples_data.to_csv(output_csv, index=False)

if __name__ == "__main__":
    input_files = list(snakemake.input.filtered_data)
    output_csv = snakemake.output.csv
    sample_name = snakemake.params.sample_name
    sys.stderr = open(snakemake.log[0], "w")  
    export_genera_abundance(input_files, sample_name, output_csv)
