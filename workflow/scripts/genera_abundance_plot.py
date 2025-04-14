import pandas as pd
import plotly.express as px
import plotly.io as pio

def create_bubble_plots(input_csv, abundance_threshold, output_html):
    df = pd.read_csv(input_csv, header=0, sep=',')

    # Drop rows with missing or zero genus counts
    df = df[df['genus_count'] > 0]

    # If empty after filtering, return an empty plot
    if df.empty or df['AMR Gene Family'].isnull().all():
        print("No valid data to plot — generating an empty dummy plot.")

        fig = px.scatter(
            title="Empty Bubble Plot — No Data Available"
        )
        fig.update_layout(
            xaxis_title='Sample - AMR Gene Family',
            yaxis_title='Genus',
            plot_bgcolor='lightgrey'
        )
        pio.write_html(fig, file=output_html)
        return

    # Normal flow
    total_counts_per_abr = df.groupby('AMR Gene Family')['genus_count'].sum().reset_index()
    top_abr = total_counts_per_abr.sort_values(by='genus_count', ascending=False)['AMR Gene Family'].iloc[0]
    top_abr_data = df[df['AMR Gene Family'] == top_abr]
    top_abr_data = top_abr_data[top_abr_data['relative_genus_count'] > float(abundance_threshold)]

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

    fig.update_layout(
        title=f'Bubble Plot of Relative Genera Abundance per Sample\nTop AMR Gene Family: {top_abr} ({total_counts_per_abr["genus_count"].sum()} reads)',
        xaxis_title='Sample - AMR Gene Family',
        yaxis_title='Genus',
        coloraxis_colorbar=dict(title="Total Genus Count"),
        plot_bgcolor='lightgrey',
        yaxis=dict(categoryorder="category descending"),
        xaxis=dict(categoryorder="category ascending"),
    )
    pio.write_html(fig, file=output_html)

if __name__ == "__main__":
    input_files = snakemake.input.abundance_data
    output_csv = snakemake.output[0]
    abundance_threshold = snakemake.params.abundance_filter
    sys.stderr = open(snakemake.log[0], "w")  
    create_bubble_plots(input_files, abundance_threshold, output_csv)
