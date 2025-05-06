import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

main1 = "Diamond and Usearch hits"
main2 = "Hits after filtration"
dia1 = "Diamond hits < similarity threshold"
dia2 = "Diamond hits ≠ max identity for query ID"
usea1 = "Usearch hits < similarity threshold"
usea2 = "Usearch hits ≠ max identity for query ID"
comb = "No overlap for hits in both databases"


def categorize_state(state):
    return {"integration_output": main1, "filtration_output": main2}.get(state)


def categorize_filter_state(state):
    return {
        "filtered_min_similarity_ABR": dia1,
        "filtered_max_identity_ABR": dia2,
        "filtered_min_similarity_16S": usea1,
        "filtered_max_identity_16S": usea2,
        "filtered_query_id_mismatch": comb,
    }.get(state)


def plot_summary(input_path, output_path):
    df = pd.read_csv(input_path)
    df["category"] = df["state"].apply(categorize_state)
    df_main = df.dropna(subset=["category"])
    main_summary = (
        df_main.groupby(["sample", "category"])["total_count"].sum().unstack().fillna(0)
    )

    df["filter_reason"] = df["state"].apply(categorize_filter_state)
    df_overlay = df.dropna(subset=["filter_reason"])
    overlay_summary = (
        df_overlay.groupby(["sample", "filter_reason"])["total_count"]
        .sum()
        .unstack()
        .fillna(0)
    )

    # Plot setup
    samples = main_summary.index
    x = np.arange(len(samples))
    bar_width = 0.25
    overlay_width = 0.125  # thinner width for filter stack

    fig, ax = plt.subplots(figsize=(12, 7))

    # Colors
    main_colors = {main1: "#fc8d62", main2: "#8da0cb"}
    filter_colors = {
        dia1: "royalblue",
        dia2: "purple",
        usea1: "#a6d854",
        usea2: "#66c2a5",
        comb: "#ffd92f",
    }

    # Plot main bars
    ax.bar(
        x - bar_width / 2,
        main_summary[main1],
        bar_width,
        label=main1,
        color=main_colors[main1],
    )
    ax.bar(
        x + bar_width / 2,
        main_summary[main2],
        bar_width,
        label=main2,
        color=main_colors[main2],
    )

    # Plot filter stack thinner, centered on same x as main2
    bottom = main_summary[main2].values.copy()
    for cat in [dia1, dia2, usea1, usea2, comb]:
        ax.bar(
            x + bar_width / 2.1,
            overlay_summary[cat],
            overlay_width,
            bottom=bottom,
            label=cat,
            color=filter_colors[cat],
        )
        bottom += overlay_summary[cat]

    # Styling
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45)
    ax.set_ylabel("Similarity search hit count")
    ax.set_xlabel("Sample")
    ax.set_title(
        "Similarity Search Processing with Rejection Breakdown on Filtration Hits"
    )

    # Legends
    handles, labels = ax.get_legend_handles_labels()
    main_labels = [main1, main2]
    filter_labels = [dia1, dia2, usea1, usea2, comb]

    main_handles = [handles[labels.index(l)] for l in main_labels]
    filter_handles = [handles[labels.index(l)] for l in filter_labels]

    legend1 = ax.legend(
        main_handles,
        main_labels,
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        title="Hit Process",
    )
    legend2 = ax.legend(
        filter_handles,
        filter_labels,
        loc="upper left",
        bbox_to_anchor=(1.02, 0.55),
        title="Filtering reasons",
    )
    ax.add_artist(legend1)

    plt.tight_layout()
    plt.savefig(output_path)


if __name__ == "__main__":
    input_path = snakemake.input.overview_table
    output_path = snakemake.output[0]
    plot_summary(input_path, output_path)
