import pandas as pd
from collections import defaultdict


def merge_overview_tables(input_paths, sample_name, output_path):
    # Dictionary of dictionaries: sample -> state -> total count
    state_counts = defaultdict(int)

    for path in input_paths:
        try:
            df = pd.read_csv(
                path, header=None, names=["state", "sample", "part", "count"]
            )
        except Exception as e:
            print(f"Could not read {path}: {e}")
            continue

        for _, row in df.iterrows():
            state = row["state"]
            count = int(row["count"])
            state_counts[state] += count

    # Write summary grouped by sample
    summary_df = pd.DataFrame(
        [
            {"sample": sample_name, "state": state, "total_count": total}
            for state, total in state_counts.items()
        ]
    )

    # Save as HTML
    summary_df.to_html(output_path, index=False)


if __name__ == "__main__":
    input_paths = snakemake.input.overview_tables
    sample_name = snakemake.params.sample_name
    output_path = snakemake.output[0]
    merge_overview_tables(input_paths, sample_name, output_path)
