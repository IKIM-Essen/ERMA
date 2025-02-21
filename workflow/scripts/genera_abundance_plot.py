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
    significant_AMRs = total_counts_per_abr[total_counts_per_abr / total_counts_per_abr.sum() > 0.2].index.tolist()

    # Create a subplot figure (no shared x/y axes or legend)
    fig = make_subplots(
        rows=len(significant_AMRs), cols=1, 
        subplot_titles=[f"{amr} ({total_counts_per_abr[amr]} reads)" for amr in significant_AMRs]
    )

    # Add traces for each AMR Gene Family
    for idx, amr in enumerate(significant_AMRs):
        amr_data = df[df['AMR Gene Family'] == amr]
        amr_data = amr_data[(amr_data['relative_genus_count'] > float(abundance_threshold)) & (amr_data['total_genus_count'] > 100)]
        amr_data = amr_data.sort_values(by=["genus", "sample"],ascending=[False,True])

        scatter = px.scatter(
            amr_data, x="sample", y="genus", 
            size="relative_genus_count", color="total_genus_count",
            hover_name="genus",
            hover_data={"genus_count": True, "relative_genus_count": True, "total_genus_count": True, "sample": False},
            size_max=20, color_continuous_scale="Greens"
        )

        # Add each scatter plot to the correct subplot
        for trace in scatter.data:
            trace.showlegend = True  # Remove legend
            fig.add_trace(trace, row=idx+1, col=1)

    # Update layout
    fig.update_layout(
        height=600 * len(significant_AMRs), 
        title="Bubble Plots of Significant AMR Gene Families Across Samples",
        plot_bgcolor='lightgrey',
        coloraxis_colorbar=dict(title="Total Genus Count")
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
