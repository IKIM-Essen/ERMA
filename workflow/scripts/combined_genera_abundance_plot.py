import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.io as pio
import sys

"""
This script creates interactive bubble plots showing the top genera per AMR Gene Family
across multiple samples based on epicPCR abundance data.

It generates one subplot per gene family with at least min_total_count hits,
displaying genus abundance as bubble size and total counts as color intensity.

Decision of which genera shown are made in the function select_genera and is whether,
the overlap of top ranking genera over all samples if the overlap is greater than
min_overlap, one of each samples top genus iteratively or every genus present.
"""


def load_filtered_data(input_csv, min_total_count=100):
    """Load CSV and filter AMR Gene Families by minimum total genus count"""
    df = pd.read_csv(input_csv, sep=",")
    return df[df["total_genus_count"] > min_total_count]


def get_top_genera_per_sample(df, top_n):
    """Return dicts of top genera per sample (set and list forms)"""
    top_sets = {}
    top_lists = {}
    for sample in df["sample"].unique():
        sample_df = df[df["sample"] == sample]
        top = (
            sample_df.sort_values(by="relative_genus_count", ascending=False)
            .head(top_n)["genus"]
            .tolist()
        )
        top_sets[sample] = set(top)
        top_lists[sample] = top
    return top_sets, top_lists


def select_genera(top_sets, top_lists, max_genera, min_overlap):
    """Select a list of genera to display using overlap or merged ranking"""
    if not top_sets:
        return []

    overlap = set.intersection(*top_sets.values())
    total_genus = sum(len(lst) for lst in top_lists.values())

    if len(overlap) >= min_overlap:
        return list(overlap)[:max_genera]
    elif total_genus > max_genera:
        combined = set(overlap)
        sample_iters = {s: iter(l) for s, l in top_lists.items()}

        while len(combined) < max_genera:
            for gen_iter in sample_iters.values():
                try:
                    while True:
                        genus = next(gen_iter)
                        if genus not in combined:
                            combined.add(genus)
                            break
                except StopIteration:
                    continue
                if len(combined) >= max_genera:
                    break
        return list(combined)
    else:
        return list({genus for sublist in top_lists.values() for genus in sublist})


def add_amr_family_subplot(
    fig, df, amr_family, col_idx, max_genera, min_overlap, top_per_sample
):
    """Filter and add a subplot for one AMR Gene Family to the main figure"""
    df_amr = df[df["AMR Gene Family"] == amr_family]
    if df_amr.empty:
        return

    top_sets, top_lists = get_top_genera_per_sample(df_amr, top_per_sample)
    selected = select_genera(top_sets, top_lists, max_genera, min_overlap)
    df_plot = df_amr[df_amr["genus"].isin(selected)]

    scatter = px.scatter(
        df_plot,
        x="sample",
        y="genus",
        size="relative_genus_count",
        color="total_genus_count",
        hover_name="genus",
        hover_data={
            "genus_count": True,
            "relative_genus_count": True,
            "total_genus_count": True,
            "sample": False,
        },
        size_max=20,
        color_continuous_scale="Greens",
    )

    for trace in scatter.data:
        fig.add_trace(trace, row=1, col=col_idx)


def create_bubble_plot_grid(df, max_genera, min_overlap, top_per_sample):
    """Create the full multi-subplot bubble chart"""
    families = df["AMR Gene Family"].unique()
    num_cols = len(families) if len(df) > 1 else 1

    fig = make_subplots(
        rows=1,
        cols=num_cols,
        subplot_titles=list(families),
        horizontal_spacing=0.2,
    )

    for idx, family in enumerate(families, start=1):
        add_amr_family_subplot(
            fig, df, family, idx, max_genera, min_overlap, top_per_sample
        )

    fig.update_layout(
        title="Bubble Plots of Top Genera for Each AMR Gene Family",
        plot_bgcolor="lightgrey",
        height=900,
        width=500 * num_cols,
        coloraxis_colorbar=dict(title="Total Genus Count"),
    )
    fig.update_yaxes(categoryorder="category descending")
    fig.update_xaxes(categoryorder="category ascending")

    return fig


def create_bubble_plots_combined(
    input_csv, output_html, max_genera=20, min_overlap=10, top_per_sample=20
):
    """Load input, pass to processing function and save plot"""
    df = load_filtered_data(input_csv)
    fig = create_bubble_plot_grid(df, max_genera, min_overlap, top_per_sample)
    pio.write_html(fig, file=output_html)


if __name__ == "__main__":
    input_csv = snakemake.input.abundance_data
    output_html = snakemake.output[0]
    sys.stderr = open(snakemake.log[0], "w")
    create_bubble_plots_combined(input_csv, output_html)
