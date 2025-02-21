import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

necessary_columns = [
    "query_id",
    "part",
    "genus",
    "AMR Gene Family",
    "perc_identity",
]

dtype_dict = {
    "query_id": "string",
    "part": "string",
    "genus": "string",
    "AMR Gene Family": "string",
    "perc_identity": "float"
}

def generate_percentage_idt_per_genus(input_files, output_file):
    all_data = []  # List to hold DataFrames from all input files
    sample_name = input_files[0].split("results/")[1].replace("/001/filtered_results.csv.gz","")
    for input_file in input_files:
        df = pd.read_csv(input_file, sep=",", usecols=necessary_columns, header=0, dtype=dtype_dict, compression='gzip')
        df = df.drop_duplicates()
        all_data.append(df)
        
    # Combine all partitions into a single DataFrame
    combined_data = pd.concat(all_data, ignore_index=True)
    
    #Filter out gene specific duplicates
    abr_data = combined_data[combined_data["part"] == "ABR"]
    sixteen_s_data = combined_data[combined_data["part"] == "16S"] 
    unique_abr_data = abr_data.drop_duplicates(['query_id', 'AMR Gene Family'])
    unique_sixteen_s_data = sixteen_s_data.drop_duplicates(['query_id', 'genus'])

    # Remerge
    merged_data = pd.merge(unique_abr_data[['query_id', 'AMR Gene Family']], unique_sixteen_s_data[['query_id', 'genus', "perc_identity"]],on='query_id', how='inner')       
    
    # Calculate genus query counts and genus order
    genus_query_counts = merged_data.groupby('genus')['query_id'].nunique().reset_index()
    genus_query_counts.columns = ['genus', 'unique_query_count']

    threshold = merged_data.shape[0]*0.005
    genus_query_counts = genus_query_counts[genus_query_counts['unique_query_count'] >= threshold]
    filtered_genus = genus_query_counts['genus'].tolist()
    merged_data = merged_data[merged_data['genus'].isin(filtered_genus)]

    genus_order = (genus_query_counts.sort_values(by='unique_query_count', ascending=False)['genus'])
    
    # Plotting
    colorboxplot = 'dodgerblue'
    colorbars = 'purple'
    
    fig, ax1 = plt.subplots(figsize=(9, 8))
    sns.boxplot(
        x='genus', y='perc_identity', data=merged_data, ax=ax1, order=genus_order, fliersize=0.0, color=colorboxplot
    )
    ax1.set_xlabel("Bacterial Genus")
    ax1.set_ylabel("Percentage Identity (boxplot)", color='royalblue')
    ax1.set_title(f"Sample {sample_name}: Bacterial genus percentage identity and read counts (threshold 0.5% of total)")
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)

    # Add a second y-axis for unique query counts
    ax2 = ax1.twinx()
    sns.barplot(
        x='genus', y='unique_query_count', data=genus_query_counts,
        ax=ax2, alpha=0.2, color=colorbars, order=genus_order
    )
    ax2.set_ylabel("Number of hits (bar)", color='violet')
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

if __name__ == "__main__":
    input_files = snakemake.input.filtered_data  # List of partitioned files
    output_file = snakemake.output[0]
    sys.stderr = open(snakemake.log[0], "w")
    generate_percentage_idt_per_genus(input_files, output_file)
