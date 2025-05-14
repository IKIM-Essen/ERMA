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
MAIN_CATEGORIES = {
    "merge output": "Diamond and Usearch hits",
    "filtration output": "Hits after filtration",
}

FILTER_REASONS = {
    "filtered min similarity ABR": "Diamond hits < similarity threshold",
    "filtered max identity ABR": "Diamond hits ≠ max identity for query ID",
    "filtered min similarity 16S": "Usearch hits < similarity threshold",
    "filtered max identity 16S": "Usearch hits ≠ max identity for query ID",
    "filtered query id mismatch": "No overlap for hits in both databases",
}


def map_main_category(state):
    return MAIN_CATEGORIES.get(state)


def map_filter_reason(state):
    return FILTER_REASONS.get(state)


def load_and_summarize_data(input_path):
    """Read the overview table and group by main and filtering categories."""
    df = pd.read_csv(input_path)

    # Assign main and filter categories
    df["category"] = df["state"].apply(map_main_category)
    df["filter_reason"] = df["state"].apply(map_filter_reason)
    df["total_count"] = df["total_count"].astype(int).abs()

    # Group main categories
    main_summary = (
        df.dropna(subset=["category"])
        .groupby(["sample", "category"])["total_count"]
        .sum()
        .unstack()
        .fillna(0)
    )

    # Group filtering reasons
    overlay_summary = (
        df.dropna(subset=["filter_reason"])
        .groupby(["sample", "filter_reason"])["total_count"]
        .sum()
        .unstack()
        .fillna(0)
    )

    return main_summary, overlay_summary


def plot_summary(main_summary, overlay_summary, output_path):
    """Generate and save a stacked bar plot showing filtering breakdown."""
    samples = main_summary.index
    x = np.arange(len(samples))
    bar_width = 0.25
    overlay_width = 0.125

    fig, ax = plt.subplots(figsize=(12, 7))

    # Define colors
    main_colors = {
        MAIN_CATEGORIES["merge output"]: "#fc8d62",
        MAIN_CATEGORIES["filtration output"]: "#8da0cb",
    }
    filter_colors = {
        FILTER_REASONS[k]: c
        for k, c in zip(
            FILTER_REASONS, ["royalblue", "purple", "#a6d854", "#66c2a5", "#ffd92f"]
        )
    }

    # Plot main bars
    ax.bar(
        x - bar_width / 2,
        main_summary[MAIN_CATEGORIES["merge output"]],
        bar_width,
        label=MAIN_CATEGORIES["merge output"],
        color=main_colors[MAIN_CATEGORIES["merge output"]],
    )
    ax.bar(
        x + bar_width / 2,
        main_summary[MAIN_CATEGORIES["filtration output"]],
        bar_width,
        label=MAIN_CATEGORIES["filtration output"],
        color=main_colors[MAIN_CATEGORIES["filtration output"]],
    )

    # Stack filter bars on top of filtration bar
    bottom = main_summary[MAIN_CATEGORIES["filtration output"]].values.copy()
    for reason in FILTER_REASONS.values():
        heights = (
            overlay_summary[reason]
            if reason in overlay_summary
            else np.zeros_like(bottom)
        )
        ax.bar(
            x + bar_width / 2.1,
            heights,
            overlay_width,
            bottom=bottom,
            label=reason,
            color=filter_colors[reason],
        )
        bottom += heights

    # Axis formatting
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45)
    ax.set_ylabel("Similarity search hit count")
    ax.set_xlabel("Sample")
    ax.set_title(
        "Similarity Search Processing with Rejection Breakdown on Filtration Hits"
    )

    # Split legend into main vs. filter categories
    handles, labels = ax.get_legend_handles_labels()
    main_labels = list(MAIN_CATEGORIES.values())
    filter_labels = list(FILTER_REASONS.values())

    legend1 = ax.legend(
        [handles[labels.index(l)] for l in main_labels],
        main_labels,
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        title="Hit Process",
    )
    legend2 = ax.legend(
        [handles[labels.index(l)] for l in filter_labels],
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
