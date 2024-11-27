import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def read_and_process_partitioned_data(partition_files,sample):
    """Read and process partitioned files for a single sample."""
    data_frames = []
    sample_name = sample

    for part_file in partition_files:
        if os.path.exists(part_file):
            fields = ['query_id','evalue', 'part']
            df = pd.read_csv(part_file, header=0, sep=',', usecols=fields)
            df = df.drop_duplicates()            
            df['sample'] = sample_name
            data_frames.append(df)
        else:
            print(f"File {part_file} does not exist.")

    if data_frames:
        return pd.concat(data_frames)
    else:
        return None

def plot_boxplots(data, output_file):
    """Plot boxplots based on the e-values for ABR and 16S parts across samples."""
    plt.figure(figsize=(15, 10))
    flierprops = dict(markerfacecolor='0.75', markersize=2, linestyle='none')
    sns.boxplot(x='sample', y='evalue', hue='part', data=data, flierprops=flierprops)
    plt.yscale('log')
    plt.title('Boxplot of e-values for ABR and 16S parts across samples -Filtered-')
    plt.xlabel('Sample')
    plt.ylabel('e-value (log scale)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def main(csv_files, sample_names, output_file):
    """Main function to process partitioned files for each sample and generate the plot."""
    all_data = []

    # Loop over each sample's partitioned CSV files
    for sample in sample_names:
        data = read_and_process_partitioned_data([file for file in csv_files if str(sample) in file],sample)
        if data is not None:
            all_data.append(data)

    if all_data:
        combined_data = pd.concat(all_data)
        print("Creating combined boxplot")
        plot_boxplots(combined_data, output_file)
        print(f"Combined boxplot has been saved to {output_file}")
    else:
        print("No data found.")

if __name__ == "__main__":
    partitioned_csv_files = snakemake.input.csv_files  # Dict of sample -> list of partition files
    output_file = snakemake.output[0]  # Path to save the output plot
    sample_name = snakemake.params.sample_name  # Minimum similarity filter
    sys.stderr = open(snakemake.log[0], "w")
    main(partitioned_csv_files, sample_name, output_file)
