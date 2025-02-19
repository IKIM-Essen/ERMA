import pandas as pd
import plotly.express as px
import sys

def create_bubble_plots(df, abundance_threshold, output1,output2):
    # Iterate over unique AMR Gene Families
    # Filter data for the current AMR Gene Family
    df = pd.read_csv(df,header=0,sep=',')
    total_counts_per_abr = df.groupby('AMR Gene Family')['genus_count'].sum()
    total_counts_per_abr.to_csv(output2)
    top_abr = total_counts_per_abr.idxmax()
    top_abr_data = df[df['AMR Gene Family'] == top_abr]
    top_abr_data = top_abr_data[top_abr_data['relative_genus_count'] > float(abundance_threshold)]
    # Create the bubble plot
    fig = px.scatter(
        top_abr_data, 
        x="sample", 
        y="genus", 
        size="relative_genus_count",
        color="total_genus_count", 
        hover_name="genus",
        hover_data={
            "genus_count": True,
            "relative_genus_count": True,
            "total_genus_count": True,
            "sample": False
        },
        size_max=20,
        color_continuous_scale="Greens"
    )
    
    # Update layout for titles and axis labels
    fig.update_layout(
        title=f'Bubble Plot of Relative Genera Abundance per Sample and top found AMR:<br> {top_abr} with {total_counts_per_abr.sum()} reads over all samples',
        xaxis_title='Sample - AMR Gene Family',
        yaxis_title='Genus',
        coloraxis_colorbar=dict(title="Total Genus Count"),
        #paper_bgcolor='grey',
        plot_bgcolor='lightgrey',
        yaxis=dict(categoryorder="category descending"),
        xaxis=dict(categoryorder="category ascending"),
    )
    fig.write_html(output1)

if __name__ == "__main__":
    input_files = snakemake.input.abundance_data
    output_html = snakemake.output[0]
    output_csv = snakemake.output[1]
    abundance_threshold = snakemake.params.abundance_filter
    sys.stderr = open(snakemake.log[0], "w")  
    create_bubble_plots(input_files, abundance_threshold, output_html,output_csv)
