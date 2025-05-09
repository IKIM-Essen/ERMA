import pandas as pd
import plotly.express as px
import sys

"""
This script takes the combined genus abundance table as input and counts
the number of hits in relation to respective AMR families.
"""


def create_bubble_plots(df, output):
    """ Group by unique AMR Gene Families and sum respective genus count """
    # Filter data for the current AMR Gene Family
    df = pd.read_csv(df, header=0, sep=",")
    total_counts_per_abr = (
        df.groupby("AMR Gene Family")["genus_count"].sum().reset_index()
    )
    total_counts_per_abr.to_html(output)


if __name__ == "__main__":
    input_files = snakemake.input.abundance_data
    output_html = snakemake.output[0]
    sys.stderr = open(snakemake.log[0], "w")
    create_bubble_plots(input_files, output_html)
