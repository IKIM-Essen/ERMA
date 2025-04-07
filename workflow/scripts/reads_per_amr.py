import pandas as pd
import plotly.express as px
import sys

def create_bubble_plots(df, output):
    # Iterate over unique AMR Gene Families
    # Filter data for the current AMR Gene Family
    df = pd.read_csv(df,header=0,sep=',')
    total_counts_per_abr = df.groupby('AMR Gene Family')['genus_count'].sum().reset_index()
    total_counts_per_abr.to_html(output)
    
if __name__ == "__main__":
    input_files = snakemake.input.abundance_data
    output_html = snakemake.output[0]
    sys.stderr = open(snakemake.log[0], "w")  
    create_bubble_plots(input_files, output_html)
