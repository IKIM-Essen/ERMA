# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


import pandas as pd
from collections import defaultdict

"""
This script merges overview table files over all samples, effectively
creating one table for all metrics between input and output.
"""


def merge_overview_tables(input_paths, output_path):
    # Dictionary of dictionaries: sample -> state -> total count
    sample_state_counts = defaultdict(lambda: defaultdict(int))

    for path in input_paths:
        try:
            df = pd.read_csv(
                path, header=None, names=["step", "sample", "part", "count"]
            )
        except Exception as e:
            print(f"Could not read {path}: {e}")
            continue

        for _, row in df.iterrows():
            sample = row["sample"]
            state = row["step"]
            count = int(row["count"])
            sample_state_counts[sample][state] += count

    # Write summary grouped by sample
    with open(output_path, "w") as out_file:
        header_written = False
        for sample, state_counts in sample_state_counts.items():
            if header_written == False:
                out_file.write("sample,step,total_count\n")
                header_written = True
            for state, total in state_counts.items():
                out_file.write(f"{sample},{state},{total}\n")


if __name__ == "__main__":
    input_paths = snakemake.input.overview_tables
    output_path = snakemake.output[0]
    merge_overview_tables(input_paths, output_path)
