# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.

rule fetch_ribosnake:
    output:
        local(touch("workflow/external/RiboSnake_reduced/.fetched")),
    run:
        import os, subprocess
        repo = config["ribosnake_repo"]
        outdir = config["ribosnake_directory"]
        
        subprocess.run(["git", "clone", repo, "-b","ERMA", outdir])


rule generate_ribosnake_metadata:
    input:
        snakefile = local("workflow/external/RiboSnake_reduced/.fetched"),
        samples = expand("data/fastq/{sample}.fastq.gz", sample=samples)
    params:
        config = config["ribosnake_config"],
        destdir = "workflow/external/RiboSnake_reduced/config/"
    output:
        config="workflow/external/RiboSnake_reduced/config/config.yaml",
        metadata="workflow/external/RiboSnake_reduced/config/pep/metadata.csv"
    run:
        import os
        import pandas as pd
        import shutil
        from datetime import date

        sample_names = [os.path.basename(x).split(".")[0] for x in input.samples]
        df = pd.DataFrame({
            "sample_name": sample_names,
            "subject": ["16S"] * len(sample_names),
            "run_date": [str(date.today())] * len(sample_names)
        })
        df.to_csv(output.metadata, index=False)
        shutil.copy(params.config,output.config)

rule run_ribosnake:
    input:
        metadata= "workflow/external/RiboSnake_reduced/config/pep/metadata.csv"
    log:
        local("logs/run_ribosnake.log")
    output:
        local(touch("results/ribosnake_done.flag")) 
    run:
        # Trigger a target inside the subworkflow
        snakemake.workflow.execute_subworkflow(
            name="ribosnake",
            targets=["results/final_output_of_ribosnake.txt"]
        )