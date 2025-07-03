import os, pathlib
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ─── Constants ─────────────────────────────────────────────────────────
RESERVED_COLOR = "rgb(217,217,217)"
AMR_MIN_FRACTION = 0.01


def get_genus_colors(all_genera):
    """Assign consistent, distinguishable colors to each genus."""
    top_colors = [
        "#D62728",  # dark red
        "#FF7F0E",  # orange
        "#8B4513",  # brown
        "#1F77B4",  # dark blue
        "#800080",  # purple
        "#7F7F7F",  # gray
        "#2CA02C",  # dark green
        "#1E90FF",  # blue
        "#BA55D3",  # medium orchid
        "#BCBD22",  # yellow-green
    ]

    fallback_palette = (
        px.colors.qualitative.Pastel
        + px.colors.qualitative.Set3
        + px.colors.qualitative.Alphabet
        + px.colors.qualitative.Light24
        + px.colors.qualitative.Bold
    )

    # Remove duplicates and reserved color from palette
    color_pool = list(dict.fromkeys(top_colors + fallback_palette))
    if RESERVED_COLOR in color_pool:
        color_pool.remove(RESERVED_COLOR)

    # Assign genera with a unique color each
    genus_list = [g for g in all_genera if g != "Others"]
    if len(genus_list) > len(color_pool):
        raise ValueError(
            f"Too many genera ({len(genus_list)}) for available color pool."
        )
    genus_colors = {g: color_pool[i] for i, g in enumerate(genus_list)}
    genus_colors["Others"] = RESERVED_COLOR
    return genus_colors


def preprocess_abundance(df, amr, min_genus_abundance, force_include, force_exclude):
    """Filter and aggregate genus abundance data for a given AMR family."""
    df_amr = df[df["AMR Gene Family"] == amr].copy()

    # Determine low-abundance or excluded genera
    low_abundance = df_amr[
        (
            (df_amr["relative_genus_count"] <= min_genus_abundance)
            & (~df_amr["genus"].isin(force_include))
        )
        | (df_amr["genus"].isin(force_exclude))
    ]
    others = (
        low_abundance.groupby(["sample", "total_count"], as_index=False)
        .agg({"relative_genus_count": "sum"})
        .assign(genus="Others")
    )
    others["sample_label"] = (
        others["sample"] + " (" + others["total_count"].astype(str) + ")"
    )

    # Remove excluded genera
    df_amr = df_amr[~df_amr["genus"].isin(force_exclude)]
    df_amr = df_amr.sort_values(
        by=["sample", "AMR Gene Family", "genus_count"], ascending=[True, False, False]
    )
    # plot high abundance or forced-includes
    df_amr_filtered = df_amr[
        (df_amr["relative_genus_count"] > min_genus_abundance)
        | (df_amr["genus"].isin(force_include))
    ]

    # Add "Others"
    df_final = pd.concat([df_amr_filtered, others], ignore_index=True)
    df_final["sample_label"] = (
        df_final["sample"] + " (" + df_final["total_count"].astype(str) + ")"
    )
    return df_final


def plot_stacked_abundance(
    observed_csv,
    output_html,
    min_genus_abundance,
    force_include=None,
    force_exclude=None,
):
    """Main function to generate a stacked bar plot of genus abundance by AMR family."""
    force_include = force_include or []
    force_exclude = force_exclude or []

    df = pd.read_csv(observed_csv)
    df = df.sort_values(["sample", "genus_count"], ascending=[True, False])

    # ─── Filter AMR families by total count ─────────────────────────────
    amr_totals = df.groupby("AMR Gene Family")["total_count"].sum()
    total_all = amr_totals.sum()
    amrs_to_plot = amr_totals[amr_totals >= total_all * AMR_MIN_FRACTION].index.tolist()

    if not amrs_to_plot:
        print("No AMR Gene Families meet the abundance threshold.")
        return

    df = df[df["AMR Gene Family"].isin(amrs_to_plot)]
    amrs = sorted(df["AMR Gene Family"].unique())
    samples = df["sample"].nunique()

    # ─── Set up subplots ────────────────────────────────────────────────
    fig = make_subplots(
        rows=len(amrs),
        cols=1,
        subplot_titles=amrs,
        shared_xaxes=True,
        vertical_spacing=0.2,
    )

    for i, amr in enumerate(amrs, start=1):
        df_amr = preprocess_abundance(
            df, amr, min_genus_abundance, force_include, force_exclude
        )
        genus_colors = get_genus_colors(df_amr["genus"].unique())

        genera = df_amr["genus"].unique()
        for genus in genera:
            genus_data = df_amr[df_amr["genus"] == genus]
            fig.add_trace(
                go.Bar(
                    x=genus_data["sample_label"],
                    y=genus_data["relative_genus_count"],
                    name=genus,
                    marker_color=genus_colors[genus],
                    showlegend=True,
                ),
                row=i,
                col=1,
            )

    # ─── Layout ────────────────────────────────────────────────────────
    fig.update_layout(
        barmode="stack",
        title="Relative Genus Abundance per AMR Gene Family",
        height=800 * len(amrs),
        width=1000 * np.log10(samples) if samples > 2 else 500,
        plot_bgcolor="white",
        yaxis=dict(tickformat=".0%"),
        legend_title="Genus",
    )

    fig.update_xaxes(tickangle=45)
    fig.update_yaxes(title_text="Relative Abundance")

    # Save and show
    fig.write_html(output_html)


if __name__ == "__main__":
    input_csv = snakemake.input.abundance_data
    output_html = snakemake.output[0]
    min_abundance = snakemake.params[0]
    sys.stderr = open(snakemake.log[0], "w")
    plot_stacked_abundance(input_csv, output_html, float(min_abundance))
