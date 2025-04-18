import pandas as pd

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

def combine_blast_data(input_file, sample_name):
    # Load and combine data from all parts for the given sample
    df = pd.read_csv(input_file, sep=",", usecols=necessary_columns, header=0, compression='gzip')
    
    # Process combined data to get genus counts and relative values
    genus_counts = process_combined_data(df, sample_name)
    return genus_counts

def export_genera_abundance(input_file, sample_name, output_html):
    sample_data = combine_blast_data(input_file, sample_name)
    sample_data = sample_data.sort_values(by=["genus_count"], ascending=False)
    sample_data.to_html(output_html, index=False)

if __name__ == "__main__":
    input_file = snakemake.input.filtered_data
    output_html = snakemake.output[0]
    sample_name = snakemake.params.sample_name
    sys.stderr = open(snakemake.log[0], "w")  
    export_genera_abundance(input_file, sample_name, output_html)
