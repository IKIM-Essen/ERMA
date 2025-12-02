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
        font-query: sans-serif;
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
    <tr><th>Sample</th><th>Uniref query</th><th>Genus</th><th>Fusion Read Count</th><th>Relative</th></tr>
</thead>
<tbody>
"""


def plot_abundance_data(input_file, html, output_html):
    df = pd.read_csv(input_file)
    grouped_sample = df.groupby("sample")
    for sample, df_sample in grouped_sample:
        grouped_query = df_sample.groupby("Uniref query")
        for query, df_query in grouped_query:
            query_rowspan = len(df_query)
            query_first_row = True
            total_reads = df_query["genus_count"].sum()
            target_line = (
                f"{query}<br><span style='font-size:0.85em'>"
                f"Total Fusion Reads: {total_reads}</span>"
            )
            for i, row in df_query.iterrows():
                html += "<tr>"
                # Sample column appears ONCE per query
                if query_first_row:
                    html += f'<td rowspan="{query_rowspan}">{sample}</td>'
                    html += f'<td rowspan="{query_rowspan}">{target_line}</td>'
                    query_first_row = False
                # Per-row data
                html += (
                    f"<td>{row['genus']}</td>"
                    f"<td>{row['genus_count']}</td>"
                    f"<td>{row['relative_genus_count']}</td>"
                )
                html += "</tr>"

    html += """
    </tbody>
    </table>
    </body>
    </html>
    """

    with open(output_html, "w") as f:
        f.write(html)


if __name__ == "__main__":
    input_file = snakemake.input.abundance_data
    output_html = snakemake.output[0]
    sys.stderr = open(snakemake.log[0], "w")
    plot_abundance_data(input_file, html, output_html)
