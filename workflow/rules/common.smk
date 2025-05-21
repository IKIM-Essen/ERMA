# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


import os


# Prevents wildcard confusion by preventing dots or slashes in sample and keeps part numerical
wildcard_constraints:
    sample="[a-zA-Z0-9_]+",
    part="\\d+",


def get_base_dir():
    return config["base_dir"]


def get_card_db_dir():
    return os.path.join(get_base_dir(), "data", "card_db")


def get_silva_db_dir():
    return os.path.join(get_base_dir(), "data", "silva_db")


def get_numpart_list():
    return [f"{i:03d}" for i in range(1, config["num_parts"] + 1)]
