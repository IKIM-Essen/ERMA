import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px


def preprocess_data(df, amr, min_total_fraction, force_include, force_exclude):
    """Prepare observed data for a specific AMR gene family."""
    df_obs = df[df["AMR Gene Family"] == amr].copy()
    threshold = min_total_fraction

    # Identify "Others" before filtering
    low_abundance = df_obs[
        ((df_obs["relative_genus_count"] <= threshold) & (~df_obs["genus"].isin(force_include))) |
        (df_obs["genus"].isin(force_exclude))
    ]
    others = (
        low_abundance.groupby(['sample', 'total_count'], as_index=False)
        .agg({"relative_genus_count": "sum"})
    )
    others["genus"] = "Others"
    others["sample_label"] = others["sample"] + " (" + others["total_count"].astype(str) + ")"

    # Remove excluded genera
    df_obs = df_obs[~df_obs["genus"].isin(force_exclude)]

    # Keep high abundance or forced-includes
    df_obs_filtered = df_obs[
        (df_obs["relative_genus_count"] > threshold) | (df_obs["genus"].isin(force_include))
    ]

    # Add "Others"
    df_final = pd.concat([df_obs_filtered, others], ignore_index=True)
    df_final["sample_label"] = df_final["sample"] + " (" + df_final["total_count"].astype(str) + ")"

    return df_final


def get_genus_colors(all_genera):
    """Generate a consistent color mapping for genera."""
    color_palette = px.colors.qualitative.D3 + px.colors.qualitative.Set3
    genus_colors = {
        genus: color_palette[i % len(color_palette)] for i, genus in enumerate(sorted(all_genera))
    }
    genus_colors["Others"] = "lightgrey"
    return genus_colors


def plot_stacked_abundance(observed_csv, output_html, min_total_fraction, force_include=None, force_exclude=None):
    """Generate a multi-subplot figure showing genus abundance for each AMR gene family."""
    if force_include is None:
        force_include = []
    if force_exclude is None:
        force_exclude = []

    df = pd.read_csv(observed_csv)
    df = df.sort_values(["sample","genus_count"],ascending=[True,False])

    # ─── Filter AMR families by total count ─────────────────────────────
    amr_totals = df.groupby("AMR Gene Family")["total_count"].sum()
    total_all = amr_totals.sum()
    amrs_to_keep = amr_totals[amr_totals >= total_all * min_total_fraction].index.tolist()
    df = df[df["AMR Gene Family"].isin(amrs_to_keep)]
    amrs = sorted(df["AMR Gene Family"].unique())
    num_amrs = len(amrs)

    if num_amrs == 0:
        print("No AMR Gene Families meet the abundance threshold.")
        return

    # ─── Predefine color palette ────────────────────────────────────────
    all_genera = df["genus"].unique().tolist() + ["Others"]
    genus_colors = get_genus_colors(all_genera)

    # ─── Set up subplots ────────────────────────────────────────────────
    fig = make_subplots(
        rows=num_amrs, cols=1,
        subplot_titles=amrs,
        shared_xaxes=False,
        shared_yaxes=True,
        vertical_spacing=0.1
    )

    for i, amr in enumerate(amrs, start=1):
        df_obs = preprocess_data(df, amr, min_total_fraction, force_include, force_exclude)

        for genus in df_obs["genus"].unique():
            genus_data = df_obs[df_obs["genus"] == genus]
            fig.add_trace(
                go.Bar(
                    x=genus_data["sample_label"],
                    y=genus_data["relative_genus_count"],
                    name=genus,
                    marker_color=genus_colors[genus],
                    showlegend=True  # Show legend for each subplot
                ),
                row=i, col=1
            )

    # ─── Layout ────────────────────────────────────────────────────────
    fig.update_layout(
        barmode="stack",
        title="Relative Genus Abundance per AMR Gene Family (total number of reads)",
        height=600 * num_amrs,
        width=1000,
        plot_bgcolor="white",
        yaxis=dict(tickformat=".0%"),
        legend_title="Genus",
    )

    fig.update_xaxes(tickangle=45)
    fig.update_yaxes(title_text="Relative Abundance")

    # Save and show
    #fig.write_html(output_html)
    fig.write_html(output_html)


if __name__ == "__main__":
    input_csv = snakemake.input.abundance_data
    output_html = snakemake.output[0]
    min_abundance = snakemake.params[0]
    sys.stderr = open(snakemake.log[0], "w")
    plot_stacked_abundance(input_csv, output_html, float(min_abundance))
