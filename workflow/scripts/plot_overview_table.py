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

# Mapping step -> State
step_to_state = {
    "Number of FastQ input reads": "Input reads",
    "Diamond output hits": "Similarity search",
    "Usearch output hits": "Similarity search",
    "Merged similarity hits": "Similarity search",
    "Diamond hits < similarity threshold": "Filtration",
    "Diamond hits NOT highest percentage identity per query": "Filtration",
    "Usearch hits < similarity threshold": "Filtration",
    "Usearch hits NOT highest percentage identity per query": "Filtration",
    "Query hit in only one of two databases": "Filtration",
    "Filtered fusion reads": "Output reads",
}


# === Load and summarize the table ===
def table_to_html(input_path, output_path):
    df = pd.read_csv(input_path, header=0)
    print(df)
    df["state"] = df["step"].map(step_to_state)

    # Reorder and sort
    df = df[["sample", "state", "step", "total_count"]]
    state_order = ["Input reads", "Similarity search", "Filtration", "Output reads"]
    df["state"] = pd.Categorical(df["state"], categories=state_order, ordered=True)
    df = df.sort_values(by=["sample", "state"])

    # === HTML with rowspan for merged cells ===

    html = """
    <html>
    <head>
    <style>
        table.styled-table {
            border-collapse: collapse;
            margin: 25px 0;
            font-size: 0.95em;
            font-family: sans-serif;
            min-width: 600px;
            box-shadow: 0 0 10px rgba(0, 0, 0, 0.15);
        }
        table.styled-table thead tr {
            background-color: #009879;
            color: #ffffff;
            text-align: left;
        }
        table.styled-table th,
        table.styled-table td {
            padding: 10px 12px;
            border: 1px solid #ddd;
        }
        table.styled-table tbody tr:nth-child(even) {
            background-color: #f3f3f3;
        }
    </style>
    </head>
    <body>
    <table class="styled-table">
    <thead>
        <tr><th>Sample</th><th>State</th><th>Step</th><th>Count</th></tr>
    </thead>
    <tbody>
    """

    # Group and track rowspans
    grouped = df.groupby(["sample", "state"], observed=False)
    for (sample, state), group in grouped:
        sample_rowspan = len(df[df["sample"] == sample])
        state_rowspan = len(group)

        first_state = True
        for i, row in group.iterrows():
            html += "<tr>"
            if i == df[df["sample"] == sample].index[0]:
                html += f'<td rowspan="{sample_rowspan}">{sample}</td>'
            if first_state:
                html += f'<td rowspan="{state_rowspan}">{state}</td>'
                first_state = False
            html += f"<td>{row['step']}</td><td>{row['total_count']}</td>"
            html += "</tr>"

    html += """
    </tbody>
    </table>
    </body>
    </html>
    """
    # Write to file
    with open(output_path, "w") as f:
        f.write(html)


if __name__ == "__main__":
    input_path = snakemake.input.overview_table
    output_path = snakemake.output[0]

    table_to_html(input_path, output_path)
