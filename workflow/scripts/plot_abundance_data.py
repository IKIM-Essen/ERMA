# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


import pandas as pd
import plotly.express as px
import sys

"""
This script takes the combined genus abundance table as input and counts
the number of hits in relation to respective AMR families.
"""

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
    table.styled-table tbody tr:hover {
        background-color: #f1f1f1;
    }
</style>
</head>
<body>
<table class="styled-table">
<thead>
    <tr><th>Sample</th><th>AMR Gene Family</th><th>Genus</th><th>Fusion Read Count</th><th>Relative</th></tr>
</thead>
<tbody>
"""


def plot_abundance_data(input_file, html, output_html):
    """Group by unique AMR Gene Families and sum respective genus count"""
    # Filter data for the current AMR Gene Family
    df = pd.read_csv(input_file, sep=",", header=0)
    grouped = df.groupby(["sample", "AMR Gene Family"])
    for (sample, family), group in grouped:
        sample_rowspan = len(df[df["sample"] == sample])
        family_rowspan = len(group)
        amr = df[(df["sample"] == sample) & (df["AMR Gene Family"] == family)]
        reads_per_amr = amr["genus_count"].sum()
        amr_line = f"{family}<br><span style='font-size: 0.85em'> Total Fusion Reads: {reads_per_amr}</span>"
        first_family = True
        for i, row in group.iterrows():
            html += "<tr>"
            if i == df[df["sample"] == sample].index[0]:
                html += f'<td rowspan="{sample_rowspan}">{sample}</td>'
            if first_family:
                html += f'<td rowspan="{family_rowspan}">{amr_line}</td>'
                first_family = False
            html += f"<td>{row['genus']}</td><td>{row['genus_count']}</td><td>{row['relative_genus_count']}</td>"
            html += "</tr>"

    html += """
    </tbody>
    </table>
    </body>
    </html>
    """
    # Write to file
    with open(output_html, "w") as f:
        f.write(html)


if __name__ == "__main__":
    input_file = snakemake.input.abundance_data
    output_html = snakemake.output[0]
    sys.stderr = open(snakemake.log[0], "w")
    plot_abundance_data(input_file, html, output_html)
