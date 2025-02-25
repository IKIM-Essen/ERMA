import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots

def create_bubble_plots(df, abundance_threshold, output1,output2,output3):
    # Iterate over unique AMR Gene Families
    # Filter data for the current AMR Gene Family
    df = pd.read_csv(df,header=0,sep=',')
    total_counts_per_abr_per_sample = df.groupby(['sample','AMR Gene Family'])['genus_count'].sum()
    total_counts_per_abr_per_sample.to_csv(output3)
    total_counts_per_abr = df.groupby('AMR Gene Family')['genus_count'].sum()
    total_counts_per_abr.to_csv(output2)
    total_genus_count = total_counts_per_abr.sum()
    high_contributing_abr = total_counts_per_abr[total_counts_per_abr / total_genus_count > 0.10].index.tolist()

    # Filter and merge data for selected AMR Gene Families
    merged_data = df[df['AMR Gene Family'].isin(high_contributing_abr)]
    merged_data = merged_data[
        (merged_data['relative_genus_count'] > float(abundance_threshold)) & 
        (merged_data['total_genus_count'] > 100)
    ]
    # Create a single combined bubble plot
    fig = px.scatter(
        merged_data, 
        x="sample", 
        y="genus",
        size="relative_genus_count",
        color="AMR Gene Family",  # Different AMR Gene Families get different colors
        hover_name="genus",
        hover_data={
            "genus_count": True,
            "relative_genus_count": True,
            "total_genus_count": True,
            "sample": False
        },
        size_max=20,
        color_discrete_sequence=px.colors.qualitative.Set1  # Better color distinction
    )
    
    # Update layout for better readability
    fig.update_layout(
        title='Bubble Plot of Relative Genera Abundance per Sample for High-Contributing AMR Gene Families',
        xaxis_title='Sample - AMR Gene Family',
        yaxis_title='Genus',
        legend_title="AMR Gene Family",
        plot_bgcolor='lightgrey',
        yaxis=dict(categoryorder="category descending"),
        xaxis=dict(categoryorder="category ascending"),
        height = 20 * len(merged_data) if len(merged_data) <= 50 else 10 * len(merged_data)
    )

    # Save the final figure as an HTML file for Snakemake
    fig.write_html(output1)

if __name__ == "__main__":
    input_files = snakemake.input.abundance_data
    output_html = snakemake.output[0]
    output_csv1 = snakemake.output[1]
    output_csv2 = snakemake.output[2]
    abundance_threshold = snakemake.params.abundance_filter
    sys.stderr = open(snakemake.log[0], "w")  
    create_bubble_plots(input_files, abundance_threshold, output_html,output_csv1,output_csv2)
