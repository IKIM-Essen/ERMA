#!/usr/bin/env python3
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# ---------------------------- #
# ---------- HELPERS ---------- #
# ---------------------------- #

def collapse_by_threshold(series, threshold, label):
    """Collapse low-frequency categories below threshold into a single label."""
    freq = series.value_counts(normalize=True)
    rare = freq[freq < float(threshold)].index
    return series.replace(rare, label)


def prepare_data(uniref_abundance, threshold):
    """Read, filter, merge and preprocess data for plotting."""
    df = pd.read_csv(uniref_abundance).rename(columns={"subject_id_ABR": "Cluster ID"})

    len_dict={}
    for query in df["Uniref query"].unique():
        count = df[df["Uniref query"] == query]
        len_dict[query]=count["genus_count"].sum()
    # Collapse rare taxa
    df["genus"] = df.groupby("Uniref query", group_keys=False)["genus"].apply(
        collapse_by_threshold, threshold=threshold, label="Silva_lowfreq_bin"
    )
    df["Common taxon"] = df.groupby("Uniref query", group_keys=False)["Common taxon"].apply(
        collapse_by_threshold, threshold=threshold, label="Uniref_lowfreq_bin"
    )
    return df,len_dict

def give_dummy_plot(output_html):
    fig = go.Figure()
    fig.add_annotation(
        text="No data available for plotting.<br>(Input files contained zero matching rows.)",
        x=0.5, y=0.5,
        xref="paper", yref="paper",
        showarrow=False,
        font=dict(size=18)
    )
    fig.update_layout(
        height=300,
        width=700,
        title_text="Fusion Read Summary (No Data)"
    )
    fig.write_html(str(output_html))
    return

# ---------------------------- #
# ---------- PLOTS ------------ #
# ---------------------------- #

def plot_summary(uniref_df,len_dict,card_abundance, output_html):
    """Generate 3-row plot grid:
    Row 1: Dummy pie chart
    Row 2: Barplot of Cluster Name frequencies
    Row 3: Sankey diagram (UniRef taxon ↔ SILVA genus)
    """

    card_df = pd.read_csv(card_abundance)
    card_count = card_df["genus_count"].sum()

    queries = uniref_df["Uniref query"].unique()
    n_queries = len(queries)
    if n_queries == 0:
        give_dummy_plot(output_html)
        return    

    row_titles = [
        "Fusion read composition by database hit",
        "Cluster Name frequency distribution",
        "Stacked barplot",
        "UniRef taxon ↔ SILVA genus mapping",
    ]

    subplot_titles = [
        f"{q}:{row_titles[r]}"
        for r in range(len(row_titles))
        for q in queries
    ]

    # Three rows: Pie → Barplot → Sankey
    fig = make_subplots(
        rows=4, cols=n_queries,
        subplot_titles=subplot_titles,
        specs=[
            [{"type": "domain"}] * n_queries,   # Row 1: pie chart
            [{"type": "xy"}] * n_queries,       # Row 2: barplot
            [{"type": "xy"}] * n_queries,      # Row 3: Stacked barplot            
            [{"type": "sankey"}] * n_queries    # Row 4: sankey
        ],
        horizontal_spacing=0.1,
        vertical_spacing=0.05,
    )

    # --- Loop per UniRef query column
    for i, q in enumerate(queries, start=1):
        df_q = uniref_df[uniref_df["Uniref query"] == q].copy()

        # ------------------------------------------------------
        # ROW 1: Pie chart (two dummy values)
        # ------------------------------------------------------
        fig.add_trace(
            go.Pie(
                labels=["CARD Hits", "Uniref Hits"],
                values=[card_count, len_dict[q]],                            
                textinfo="label+percent",                
                hoverinfo="label+percent",                
                hole=0.3
            ),
            row=1, col=i
        )

        # ------------------------------------------------------
        # ROW 2: Barplot of Cluster Name frequencies
        # ------------------------------------------------------
        freq = df_q["Cluster Name"].value_counts()
        fig.add_trace(
            go.Bar(
                x=freq.values,
                y=freq.index,
                orientation="h",
                marker=dict(opacity=0.75),
                name=f"{q} frequencies"
            ),
            row=2, col=i
        )

        fig.update_yaxes(
            automargin=True,
            row=2, col=i
        )

        # ------------------------------------------------------
        # ROW 3: Genus abundance
        # ------------------------------------------------------
        genus_counts = (
            df_q["genus"].value_counts(normalize=True)
            .reset_index()
            .rename(columns={"proportion": "relative_abundance"})
        )
        for _, row in genus_counts.iterrows():
            fig.add_trace(
                go.Bar(
                    x=[q],
                    y=[row["relative_abundance"]],
                    name=row["genus"],
                    width=0.2,
                ),
                row=3, col=i
            )

        fig.update_layout(
            barmode="stack")

        # ------------------------------------------------------
        # ROW 4: Sankey diagram
        # ------------------------------------------------------
        genus_labels = list(df_q["genus"].unique())
        taxon_labels = list(df_q["Common taxon"].unique())
        all_labels = genus_labels + taxon_labels
        idx = {k: v for v, k in enumerate(all_labels)}

        sankey_links = {
            "source": [idx[t] for t in df_q["Common taxon"]],
            "target": [idx[g] for g in df_q["genus"]],
            "value": [1] * len(df_q)
        }

        sankey = go.Sankey(
            node=dict(label=all_labels, pad=10, thickness=12),
            link=sankey_links
        )
        fig.add_trace(sankey, row=4, col=i)

    # Layout ----------------------------------------------------
    fig.update_layout(
        height=1700,
        width=450 * n_queries,
        title_text=f"Fusion Read Summary of non-CARD targets: {', '.join(queries)}",
        showlegend=False
    )

    for ann in fig['layout']['annotations']:
        ann['yshift'] = 10          # move title upward (increase for more space)
        ann['font'] = dict(size=14) 

    # Save output
    fig.write_html(str(output_html))


# ---------------------------- #
# ---------- MAIN ------------- #
# ---------------------------- #

def main():
    """Entry point for Snakemake execution."""
    uniref_abundance = snakemake.input.uniref_abundance
    card_abundance = snakemake.input.card_abundance
    output_html = snakemake.output
    threshold = snakemake.params.low_freq_threshold 

    uniref_df,len_dict = prepare_data(uniref_abundance, threshold)
    plot_summary(uniref_df,len_dict,card_abundance, output_html)


if __name__ == "__main__":
    main()
