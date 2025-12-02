import requests
import pandas as pd
import time
import io

card_fasta = snakemake.input.card_fasta
output_fasta = snakemake.output.merged_fasta
output_tsv = snakemake.output.info_tsv
targets = snakemake.params.targets
cluster = snakemake.params.cluster
size = int(snakemake.params.size)
log = snakemake.log[0]

base_url = "https://rest.uniprot.org/uniref/search"
headers = {"accept": "text/plain"}
identity_map = {"100": "1.0", "90": "0.9", "50": "0.5"}


def fetch_all(base_url, params, headers, max_entries):
    """Fetch pages from UniProt REST API with an optional total entry limit."""
    all_data = []
    total_entries = 0

    while True:
        resp = requests.get(base_url, headers=headers, params=params)
        resp.raise_for_status()
        text = resp.text

        # Count FASTA entries
        n_entries = text.count(">")
        total_entries += n_entries
        all_data.append(text)

        # Stop if limit reached
        if max_entries and total_entries >= max_entries:
            print(f"Reached max_entries={max_entries}, stopping.")
            break

        # Check for next page
        next_link = resp.links.get("next", {}).get("url")
        if not next_link:
            break
        base_url = next_link
        params = {}
        time.sleep(0.5)

    return "".join(all_data)


def df_to_fasta(df):
    """Convert DataFrame to FASTA formatted string."""
    fasta_lines = []
    for _, row in df.iterrows():
        header = f">{row['Cluster ID']} {row['Common taxon']}".strip()
        seq = row["Reference sequence"].replace(" ", "")
        fasta_lines.append(f"{header}\n{seq}")
    return "\n".join(fasta_lines)


df_final = []
for t in targets:
    query = f"{t} AND identity:{identity_map[cluster]}"
    params = {
        "query": query,
        "fields": ["id", "name", "common_taxon", "identity", "sequence"],
        "format": "tsv",
    }
    fetch = fetch_all(base_url, params, headers, size)
    df = pd.read_csv(io.StringIO(fetch), header=0, sep="\t")
    df["Uniref query"] = t
    df["db"] = "Uniref"
    df_final.append(df)
if df_final:
    result = (
        pd.concat(df_final, ignore_index=True).drop_duplicates().reset_index(drop=True)
    )
    result["Cluster Name"] = result["Cluster Name"].str.replace("Cluster:", "")
    result = result.drop(
        result[result["Common taxon"] == "Common taxon"].index
    )  # remove double headers
    result.to_csv(output_tsv, index=False)
    fasta_text = df_to_fasta(result)
    with open(card_fasta, "r") as cf, open(output_fasta, "w") as out:
        out.write(cf.read().strip() + "\n" + fasta_text)

    print(f"Final combined file written: {output_fasta}")
else:
    print("No data fetched.")
