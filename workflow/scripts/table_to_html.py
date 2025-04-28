# ../scripts/table_to_html.py

import pandas as pd

def main():
    input_file = snakemake.input[1]  # overview_tables
    output_file = snakemake.output[0]

    # Read the input table with comma separator
    df = pd.read_csv(input_file, sep=",")

    # Convert to HTML and save
    df.to_html(output_file, index=False, escape=False)

if __name__ == "__main__":
    main()
