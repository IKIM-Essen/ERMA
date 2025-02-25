import os

def get_base_dir():
    return config["base_dir"]

def get_fastq_dir():
    return os.path.join(get_base_dir(), "data", "fastq")

def get_blast_db_dir():
    return os.path.join(get_base_dir(), "data", "blast_db")

def get_output_dir():
    return os.path.join(get_base_dir(), "results")

def get_card_db_dir():
    return os.path.join(get_base_dir(), "data", "card_db")

def get_silva_db_dir():
    return os.path.join(get_base_dir(), "data", "silva_db")

def get_numpart_list():
    return ["{:03d}".format(i) for i in range(1, config["num_parts"] + 1)]