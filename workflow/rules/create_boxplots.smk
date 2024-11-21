reads = [os.path.basename(f).replace(".fastq.gz", "") for f in glob.glob(os.path.join(config["base_dir"], "data", "fastq", "*.fastq.gz"))]
# Retrieve sample names. For Illumina this is done by splitting at _R1_/_R2_ ...
if config["seq_tech"] == "Illumina":
    samples = list(set(re.split(r'_R\d_', r)[0] for r in reads))
elif config["seq_tech"] == "ONT":
    samples = reads
else:
    raise ValueError("Invalid sequencing technology specified. Check config file and README.")

rule generate_percidt_genus:
    input:
        filtered_data = expand("{{base_dir}}/results/{{sample}}/{part}/filtered_results.csv", 
                                part=get_numpart_list())
    output:
        report(
            "{base_dir}/results/{sample}/genus_idt_per_genus_plot.png",
            caption = "../../report/genus_top_hits.rst",
            category="2. Genus percentage Identity",
        )
    log:
        "{base_dir}/logs/generate_percidt_genus/{sample}.log"            
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/percidt_per_genus.py"

rule plot_alignment_length_boxplot:
    input:
        csv_files = expand("{{base_dir}}/results/{sample}/{part}/filtered_results.csv",sample=samples,part=get_numpart_list())
    output:
        report(
            "{base_dir}/results/boxplots/combined_allength_boxplot.png",
            caption = "../../report/genus_top_hits.rst",
            category="3. General data",        
        )        
    params:
        sample_name = samples,
    log:
        "{base_dir}/logs/plot_alignment_length_boxplot/combined.log"               
    conda:
        "../envs/python.yaml"     
    script:
        "../scripts/align_lengths_boxplots.py"

rule plot_percentage_identity_boxplot:
    input:
        csv_files = expand("{{base_dir}}/results/{sample}/{part}/filtered_results.csv",sample=samples,part=get_numpart_list())
    output:
        report(
            "{base_dir}/results/boxplots/combined_percidt_boxplot.png",
            caption = "../../report/genus_top_hits.rst",
            category="3. General data",        
        )
    params:
        sample_name = samples,
    log:
        "{base_dir}/logs/plot_percentage_identity_boxplot/combined.log"    
    conda:
        "../envs/python.yaml"     
    script:
        "../scripts/percidt_boxplots.py"

rule plot_evalue_boxplot:
    input:
        csv_files = expand("{{base_dir}}/results/{sample}/{part}/filtered_results.csv",sample=samples,part=get_numpart_list())
    output:
        report(
            "{base_dir}/results/boxplots/combined_evalue_boxplot.png",
            caption = "../../report/genus_top_hits.rst",
            category="3. General data"
        )
    params:
        sample_name = samples,
    log:
        "{base_dir}/logs/plot_evalue_boxplot/combined.log"               
    conda:
        "../envs/python.yaml"     
    script:
        "../scripts/evalue_boxplots.py"