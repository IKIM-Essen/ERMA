# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

"""
This script generates a summary bar plot illustrating the impact of filtration steps
on BLAST similarity hits (from Diamond and Usearch) across multiple samples.

It reads a summary CSV with hit counts and states, 
groups the data by sample and category, and visualizes both retained and
filtered hits with color-coded stacked bars.
"""

# === Constants for category mapping ===
MAIN_CATEGORIES = [
    "Number of FastQ input reads",
    "Merged similarity hits",
    "Filtered fusion reads",
]

FILTER_REASONS = {
    "Diamond hits < similarity threshold": "royalblue",
    "Diamond hits NOT highest percentage identity per query": "purple",
    "Usearch hits < similarity threshold": "#a6d854",
    "Usearch hits NOT highest percentage identity per query": "#66c2a5",
    "Query hit in only one of two databases": "#ffd92f",
}

MAIN_COLOR_MAP = {
    "Number of FastQ input reads": "seagreen",
    "Merged similarity hits": "#fc8d62",
    "Filtered fusion reads": "#8da0cb",
}

# === Load and summarize the table ===
def load_and_summarize_data(path):
    df = pd.read_csv(path, header=0)
    df["total_count"] = df["total_count"].astype(int).abs()

    main_df = df[df["step"].isin(MAIN_CATEGORIES)].pivot(index="sample", columns="step", values="total_count").fillna(0)
    filter_df = df[df["step"].isin(FILTER_REASONS)].pivot(index="sample", columns="step", values="total_count").fillna(0)

    return main_df, filter_df

# === Plotting function ===
def plot_summary(main_df, filter_df, output_path):
    samples = main_df.index
    x = np.arange(len(samples))
    bar_width = 0.18
    overlay_width = 0.1

    fig, ax = plt.subplots(figsize=(12, 7))

    # Plot main bars with offsets
    offsets = np.linspace(-bar_width, bar_width, len(MAIN_CATEGORIES))
    for i, col in enumerate(MAIN_CATEGORIES):
        if col not in main_df.columns:
            continue
        ax.bar(
            x + offsets[i],
            main_df[col],
            bar_width,
            label=col,
            color=MAIN_COLOR_MAP.get(col, "gray"),
        )

    # Plot filter stack bars *on top* of "Filtered fusion reads"
    if "Filtered fusion reads" in main_df.columns:
        bottom = main_df["Filtered fusion reads"].values.copy()
    else:
        bottom = np.zeros_like(x)

    for reason in FILTER_REASONS:
        heights = filter_df[reason].values if reason in filter_df.columns else np.zeros_like(x)
        ax.bar(
            x + bar_width,
            heights,
            overlay_width,
            bottom=bottom,
            label=reason,
            color=FILTER_REASONS.get(reason, "gray"),
        )
        bottom += heights

    # Axis formatting
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45)
    ax.set_ylabel("Hit count")
    ax.set_xlabel("Sample")
    ax.set_title("Similarity Search Processing with Rejection Breakdown")

    # Split legend into main vs. filter
    handles, labels = ax.get_legend_handles_labels()
    main_labels = MAIN_CATEGORIES
    filter_labels = FILTER_REASONS

    legend1 = ax.legend(
        [handles[labels.index(l)] for l in main_labels if l in labels],
        main_labels,
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        title="Hit Process",
    )
    legend2 = ax.legend(
        [handles[labels.index(l)] for l in filter_labels if l in labels],
        filter_labels,
        loc="upper left",
        bbox_to_anchor=(1.02, 0.55),
        title="Filtering Reasons",
    )
    ax.add_artist(legend1)

    plt.tight_layout()
    plt.savefig(output_path)


if __name__ == "__main__":
    input_path = snakemake.input.overview_table
    output_path = snakemake.output[0]

    main_summary, overlay_summary = load_and_summarize_data(input_path)
    plot_summary(main_summary, overlay_summary, output_path)
