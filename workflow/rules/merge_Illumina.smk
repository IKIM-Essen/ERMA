rule prepare_fastqs:
    output:
        sentinel=local(temp("data/fastq/merged_done.txt")),
    conda:
        "../envs/python.yaml"    
    script:
        "../scripts/prepare_fastq.sh"