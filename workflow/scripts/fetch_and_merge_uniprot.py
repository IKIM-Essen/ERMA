import requests
from pathlib import Path
import time

card_fasta = Path(snakemake.input.card_fasta)
output_fasta = Path(snakemake.output.merged_fasta)
targets = snakemake.params.targets
cluster = snakemake.params.cluster
size = int(snakemake.params.size)
log = Path(snakemake.log[0])

base_url = "https://rest.uniprot.org/uniref/search"
headers = {"accept": "text/plain"}
identity_map = {100: "1.0", 90: "0.9", 50: "0.5"}

def fetch_all(base_url, params, headers):
    """Fetch all pages from UniProt REST API."""
    response = requests.get(base_url, headers=headers, params=params)
    response.raise_for_status()
    return response.text

with open(log,"w") as lf, open(output_fasta,"w") as out:
    lf.write(f"Extending CARD with UniRef{cluster} sequences for targets: {targets}\n")

    # Write original CARD sequences first
    with open(card_fasta, "r") as cf:
        out.write(cf.read())

    for t in targets:
        query = f"{t} AND identity:{identity_map[cluster]}"
        params = {
            "query": query,
            "format": "fasta",
            "size": size
            }
        lf.write(f"Fetching UniRef{cluster} sequences for: {t}\n")
        try:
            fasta = fetch_all(base_url,params,headers)
            out.write(fasta)
        except Exception as e:
            lf.write(f"Failed to fetch {t}: {e}\n")
        time.sleep(0.5) # UniProt’s API enforces ~3 requests/sec.
