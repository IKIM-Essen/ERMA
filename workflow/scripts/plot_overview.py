import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def categorize_output_state(state):
    return {
        "fasta_input": "Raw Fasta Input",
        "integration_output": "Diamond and Usearch hits",
        "filtration_output": "Hits after filtration"
    }.get(state)

def categorize_filter_state(state):
    return {
        "filtered_min_similarity_ABR": "ABR < similarity threshold",
        "filtered_min_similarity_16S": "16S < similarity threshold",
        "filtered_query_id_mismatch": "Query ID mismatch"
    }.get(state)

def plot_summary(input_path, output_path):
    df = pd.read_csv(input_path)
    df["category"] = df["state"].apply(categorize_output_state)
    df_main = df.dropna(subset=["category"])
    main_summary = df_main.groupby(["sample", "category"])["total_count"].sum().unstack().fillna(0)

    # Overlay data
    df["filter_reason"] = df["state"].apply(categorize_filter_state)
    df_overlay = df.dropna(subset=["filter_reason"])
    overlay_summary = df_overlay.groupby(["sample", "filter_reason"])["total_count"].sum().unstack().fillna(0)

    # Reorder columns
    main_order = ["Raw Fasta Input", "Diamond and Usearch hits", "Hits after filtration"]
    filter_order = ["ABR < similarity threshold", "16S < similarity threshold", "Query ID mismatch"]

    main_summary = main_summary[main_order]
    overlay_summary = overlay_summary[filter_order]

    # Plotting
    fig, ax = plt.subplots(figsize=(14, 7))
    samples = main_summary.index
    x = np.arange(len(samples))
    bar_width = 0.3
    overlay_width = 0.1
    space = 0.01  # space between bars

    # Colors
    main_colors = {
        "Raw Fasta Input": "seagreen",
        "Diamond and Usearch hits": "#fc8d62",
        "Hits after filtration": "#8da0cb"
    }
    filter_colors = {
        "ABR < similarity threshold": "royalblue",
        "16S < similarity threshold": "#a6d854",
        "Query ID mismatch": "#ffd92f"
    }

    # Positions for main bars
    positions = {
        "Raw Fasta Input": x - bar_width - space,
        "Diamond and Usearch hits": x,
        "Hits after filtration": x + bar_width + space,
    }
    # Position for thin filter overlay
    overlay_x = x - bar_width - overlay_width + 40*space

    # Plot main (non-stacked) bars
    for cat in main_order:
        ax.bar(positions[cat], main_summary[cat], bar_width, label=cat, color=main_colors[cat])

    # Plot stacked overlay bars
    bottom_overlay = np.zeros(len(overlay_summary))
    for cat in filter_order:
        ax.bar(overlay_x, overlay_summary[cat], overlay_width, bottom=bottom_overlay, label=cat, color=filter_colors[cat])
        bottom_overlay += overlay_summary[cat]

    # Styling
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45)
    ax.set_ylabel("Total Count")
    ax.set_xlabel("Sample")
    ax.set_title("Sample Processing Overview with Filtering Breakdown")

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    main_labels = [l for l in main_order]
    filter_labels = [l for l in filter_order]

    # Group legend entries
    main_handles = [handles[labels.index(l)] for l in main_labels]
    filter_handles = [handles[labels.index(l)] for l in filter_labels]

    first_legend = ax.legend(main_handles, main_labels, loc="upper left", bbox_to_anchor=(1.05, 1), title="Main Categories")
    second_legend = ax.legend(filter_handles, filter_labels, loc="upper left", bbox_to_anchor=(1.05, 0.5), title="Filter Reasons")
    ax.add_artist(first_legend)

    plt.tight_layout()
    plt.savefig(output_path)

if __name__ == "__main__":
    input_path = snakemake.input.overview_table
    output_path = snakemake.output[0]
    plot_summary(input_path,output_path)