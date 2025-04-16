import pandas as pd
import plotly.express as px
import plotly.io as pio

def create_bubble_plots(input_csv, output_html, max_genera=20, min_overlap=10, top_per_sample=20):
    df = pd.read_csv(input_csv, sep=',')

    # Sanity check: skip if empty
    if df.empty or 'genus_count' not in df.columns:
        print("No valid data to plot.")
        fig = px.scatter(title="No Data Available")
        pio.write_html(fig, file=output_html)
        return

    # Identify top AMR Gene Family by total genus_count
    top_abr = df.groupby('AMR Gene Family')['genus_count'].sum().idxmax()
    top_abr_data = df[df['AMR Gene Family'] == top_abr]

    # Get list of samples
    samples = top_abr_data['sample'].unique()

    # Collect top N genera per sample
    top_genera_per_sample_set = {}
    top_genera_per_sample_list = {}
    for sample in samples:
        sample_df = top_abr_data[top_abr_data['sample'] == sample]
        top_genera = sample_df.sort_values(by='relative_genus_count', ascending=False).head(top_per_sample)['genus'].tolist()
        top_genera_per_sample_set[sample] = set(top_genera)
        top_genera_per_sample_list[sample] = top_genera

    # Find overlap across all samples
    genus_overlap = set.intersection(*top_genera_per_sample_set.values())

    # Decide which genera to keep
    if len(genus_overlap) >= min_overlap:
        selected_genera = list(genus_overlap)[:max_genera]
    else:
        combined_genera = set(genus_overlap)
        sample_iters = {sample: iter(genera) for sample, genera in top_genera_per_sample_list.items()}
        # Keep adding unique genera in round-robin fashion until we hit the desired number
        while len(combined_genera) < max_genera:
            for sample, gen_iter in sample_iters.items():
                try:
                    while True:
                        genus = next(gen_iter)
                        if genus not in combined_genera:
                            combined_genera.add(genus)
                            break  # add only one new genus per sample per round
                except StopIteration:
                    continue  # no more genera left in this sample's list
                if len(combined_genera) >= max_genera:
                    break

        selected_genera = list(combined_genera)

    # Filter to selected genera
    plot_data = top_abr_data[top_abr_data['genus'].isin(selected_genera)]

    # Plot
    fig = px.scatter(
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
        color_continuous_scale="Greens"
    )

    fig.update_layout(
        title=f'Bubble Plot of Top Genera for AMR Gene Family: {top_abr}',
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
    sys.stderr = open(snakemake.log[0], "w")  
    create_bubble_plots(input_files, output_csv)
