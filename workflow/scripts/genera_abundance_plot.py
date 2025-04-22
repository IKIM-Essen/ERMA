import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.io as pio

def create_bubble_plots_combined(input_csv, output_html, max_genera=20, min_overlap=10, top_per_sample=20):
    df = pd.read_csv(input_csv, sep=',')

    # Filter AMR gene families with sufficient total genus count
    df = df[df["total_genus_count"] > 100]
    amr_gene_families = df['AMR Gene Family'].unique()
    num_families = len(amr_gene_families)

    # Create subplots
    fig = make_subplots(
        rows=1, cols=num_families,
        subplot_titles=[f"{amr}" for amr in amr_gene_families],
        horizontal_spacing=0.2
    )

    for i, amr in enumerate(amr_gene_families, start=1):
        abr_data = df[df['AMR Gene Family'] == amr]

        if abr_data.empty:
            continue

        samples = abr_data['sample'].unique()

        # Collect top N genera per sample
        top_genera_per_sample_set = {}
        top_genera_per_sample_list = {}
        for sample in samples:
            sample_df = abr_data[abr_data['sample'] == sample]
            top_genera = sample_df.sort_values(by='relative_genus_count', ascending=False).head(top_per_sample)['genus'].tolist()
            top_genera_per_sample_set[sample] = set(top_genera)
            top_genera_per_sample_list[sample] = top_genera

        # Find overlap across all samples
        genus_overlap = set.intersection(*top_genera_per_sample_set.values()) if top_genera_per_sample_set else set()

        # Decide which genera to keep
        if len(genus_overlap) >= min_overlap:
            selected_genera = list(genus_overlap)[:max_genera]
        else:
            combined_genera = set(genus_overlap)
            sample_iters = {sample: iter(genera) for sample, genera in top_genera_per_sample_list.items()}
            while len(combined_genera) < max_genera:
                for sample, gen_iter in sample_iters.items():
                    try:
                        while True:
                            genus = next(gen_iter)
                            if genus not in combined_genera:
                                combined_genera.add(genus)
                                break
                    except StopIteration:
                        continue
                    if len(combined_genera) >= max_genera:
                        break
            selected_genera = list(combined_genera)

        plot_data = abr_data[abr_data['genus'].isin(selected_genera)]

        scatter = px.scatter(
            plot_data, 
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
            color_continuous_scale="Greens",
        )

        for trace in scatter.data:
            fig.add_trace(trace, row=1, col=i)

    fig.update_layout(
        title='Bubble Plots of Top Genera for Each AMR Gene Family',
        plot_bgcolor='lightgrey',
        height=900,
        width=500 * num_families,
        coloraxis_colorbar=dict(title="Total Genus Count")
    )
    fig.update_yaxes(categoryorder="category descending")
    fig.update_xaxes(categoryorder="category ascending")
    pio.write_html(fig, file=output_html)

if __name__ == "__main__":
    input_files = snakemake.input.abundance_data
    output_csv = snakemake.output[0]
    sys.stderr = open(snakemake.log[0], "w")  
    create_bubble_plots_combined(input_files, output_csv)
