import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

"""
This script visualizes the genus-wise percentage identity and the number of found hits per sample.
It takes the csv files containing ABR and 16S filtered data as input. The one plot output is:

- A boxplot showing percentage identity per genus (left y-axis)
- A translucent barplot overlaying the count of unique query IDs (right y-axis)
"""

# Required columns and dtypes for input data
NECESSARY_COLUMNS = ["query_id", "part", "genus", "AMR Gene Family", "perc_identity"]
DTYPE_DICT = {
    "query_id": "string",
    "part": "string",
    "genus": "string",
    "AMR Gene Family": "string",
    "perc_identity": "float",
}


def load_data(input_files):
    """Load and concatenate all filtered CSV files."""
    dataframes = [
        pd.read_csv(
            f,
            sep=",",
            usecols=NECESSARY_COLUMNS,
            header=0,
            dtype=DTYPE_DICT,
            compression="gzip",
        )
        for f in input_files
    ]
    return pd.concat(dataframes, ignore_index=True)


def get_top_n_genera(df, n=20):
    """Return filtered data for top-N genera by unique query_id count."""
    genus_counts = df.groupby("genus")["query_id"].nunique().reset_index()
    genus_counts.columns = ["genus", "unique_query_count"]
    top_genera = genus_counts.nlargest(n, "unique_query_count")

    # Filter both combined data and counts
    top_data = df[df["genus"].isin(top_genera["genus"])]
    top_counts = genus_counts[genus_counts["genus"].isin(top_genera["genus"])]

    # Define genus order for plotting
    genus_order = top_genera.sort_values("unique_query_count", ascending=False)["genus"]

    return top_data, top_counts, genus_order


def plot_identity_boxplot(data, genus_counts, genus_order, output_file):
    """Create and save the composite plot showing identity and query abundance."""
    fig, ax1 = plt.subplots(figsize=(15, 8))

    # Boxplot: Percent identity
    sns.boxplot(
        x="genus",
        y="perc_identity",
        data=data,
        ax=ax1,
        order=genus_order,
        fliersize=0.0,
        color="dodgerblue",
    )
    ax1.set_xlabel("Bacterial Genus")
    ax1.set_ylabel("Percentage Identity (boxplot)", color="royalblue")
    ax1.set_title("Boxplot of Percentage Identity and Read Counts for Each Bacterial Genus")
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)

    # Barplot: Unique query count on secondary y-axis
    ax2 = ax1.twinx()
    sns.barplot(
        x="genus",
        y="unique_query_count",
        data=genus_counts,
        ax=ax2,
        alpha=0.2,
        color="purple",
        order=genus_order,
    )
    ax2.set_ylabel("Number of Hits (bar)", color="violet")

    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


def generate_percentage_idt_per_genus(input_files, output_file):
    """Main routine: Load data, process top genera, and create plot."""
    combined_data = load_data(input_files)
    top_data, genus_counts, genus_order = get_top_n_genera(combined_data)
    plot_identity_boxplot(top_data, genus_counts, genus_order, output_file)


if __name__ == "__main__":
    input_files = sorted(snakemake.input.filtered_data)  # List of partitioned files
    output_file = snakemake.output[0]

    # Redirect stderr to log file for Snakemake logging
    sys.stderr = open(snakemake.log[0], "w")

    generate_percentage_idt_per_genus(input_files, output_file)
