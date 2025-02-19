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

    for input_file in input_files:
        df = pd.read_csv(input_file, sep=",", usecols=necessary_columns, header=0, dtype=dtype_dict, compression='gzip')
        all_data.append(df)
        
    # Combine all partitions into a single DataFrame
    combined_data = pd.concat(all_data)
    
    # Calculate genus query counts and genus order
    genus_query_counts = combined_data.groupby('genus')['query_id'].nunique().reset_index()
    genus_query_counts.columns = ['genus', 'unique_query_count']
    combined_data = pd.merge(combined_data, genus_query_counts, on='genus')
    genus_order = genus_query_counts.sort_values(by='unique_query_count', ascending=False)['genus']
    
    # Plotting
    fig, ax1 = plt.subplots(figsize=(15, 8))
    sns.boxplot(x='genus', y='perc_identity', data=combined_data, ax=ax1, order=genus_order, fliersize=0.0)
    ax1.set_xlabel("Bacterial Genus")
    ax1.set_ylabel("Percentage Identity")
    ax1.set_title("Boxplot of Percentage Identity and Read Counts for Each Bacterial genus")
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)

    # Add a second y-axis for unique query counts
    ax2 = ax1.twinx()
    sns.barplot(x='genus', y='unique_query_count', data=genus_query_counts, ax=ax2, alpha=0.3, color='blue', order=genus_order)
    ax2.set_ylabel("Number of Unique Reads (query_id)", color='blue')
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

if __name__ == "__main__":
    input_files = snakemake.input.filtered_data  # List of partitioned files
    output_file = snakemake.output[0]
    sys.stderr = open(snakemake.log[0], "w")
    generate_percentage_idt_per_genus(input_files, output_file)
