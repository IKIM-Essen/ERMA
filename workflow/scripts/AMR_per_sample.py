import pandas as pd
import matplotlib.pyplot as plt

def main(abundance_data, output_file):
    """Main function to process partitioned files for each sample and generate the plot."""
    df = pd.read_csv(abundance_data,header=0,sep=',')
    df = df[df["genus_count"] >= 0.05 * df.groupby("sample")["genus_count"].transform("sum")]
    df_pivot = df.pivot_table(index="sample", columns="AMR Gene Family", values="genus_count", aggfunc="sum").fillna(0)

    # Plot the stacked bar chart
    fig, ax = plt.subplots(figsize=(12, 6))
    df_pivot.plot(kind="bar", stacked=True, colormap="tab10", ax=ax)

    # Customize plot
    ax.set_ylabel("Genus Count")
    ax.set_xlabel("Sample")
    ax.set_title("Stacked Bar Plot of Genus Counts by AMR Gene Family")
    ax.legend(title="AMR Gene Family", bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.xticks(rotation=60)
    plt.tight_layout()
    plt.savefig(output_file)

if __name__ == "__main__":
    abundance_data = snakemake.input.abundance_data
    output_png = snakemake.output[0]
    sys.stderr = open(snakemake.log[0], "w")
    main(abundance_data, output_png)
