import pandas as pd
import altair as alt
import sys

# Define the columns to be imported and their dtypes
necessary_columns = [
    "query_id",
    "part",
    "genus",
    "AMR Gene Family"
]

def generate_genus_distribution_plot(input_files, output_file,output2, sample_name):
    all_data = []  # List to hold DataFrames from all input files
    j=0
    for input_file in input_files:
        df = pd.read_csv(input_file, sep=",", usecols=necessary_columns, header=0, compression='gzip')
        df = df.drop_duplicates()
        # Split and merge ABR and 16S data by query_id
        abr = df[df["part"] == "ABR"]
        sixteen_s = df[df["part"] == "16S"]
        df_merged = pd.merge(abr, sixteen_s, on='query_id', suffixes=('_abr', '_16S'))
        
        # Assign genus based on 16S path for each merged record
        j+=df_merged.shape[0]
        all_data.append(df_merged)  # Append the DataFrame to the list
        #print("all_data",sample_name,j)
    # Concatenate all DataFrames into one
    combined_df = pd.concat(all_data, ignore_index=True)
    print("combined_df",sample_name,combined_df.shape[0])
    
    # Group by AMR gene family and genus across all files
    genus_distribution = combined_df.groupby(['AMR Gene Family_abr', 'genus_16S']).size().reset_index(name='count')
    genus_distribution.to_csv(output2)
    # Create the Altair bar chart
    chart = alt.Chart(genus_distribution).mark_bar().encode(
        x=alt.X('AMR Gene Family_abr', title='AMR Gene Family'),
        y=alt.Y('count', title='Count'),
        color=alt.Color('genus_16S', title='Genus'),
        tooltip=['AMR Gene Family_abr', 'genus_16S', 'count']
    ).properties(
        title=f'Genus Distribution for {sample_name}',
        width=800,
        height=400
    )
    
    # Save the plot as an HTML file
    chart.save(output_file)

input_files = snakemake.input.filtered_data
output_file = snakemake.output[0]
output2 = snakemake.output[1]
sample_name = snakemake.params.sample_name
sys.stderr = open(snakemake.log[0], "w")
print(output_file,sample_name)
generate_genus_distribution_plot(input_files, output_file,output2, sample_name)
