# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


import gzip
import re
import sys
import os


def validate_samples(samples, fastq_dir):
    sample_pattern = re.compile(r"^[a-zA-Z0-9_]+$")

    errors_found = False

    for sample in samples:
        fastq_path = os.path.join(fastq_dir, f"{sample}.fastq.gz")
        sample_errors = []

        # 1. Check if file exists
        if not os.path.isfile(fastq_path):
            sample_errors.append(f"File not found: {fastq_path}")
            errors_found = True
            print_errors(sample, sample_errors)
            continue

        # 2. Check if file is gzip-compressed
        try:
            with gzip.open(fastq_path, "rb") as f:
                f.read(1)
        except (OSError, gzip.BadGzipFile):
            sample_errors.append("File is not a valid gzip file.")
            errors_found = True
            print_errors(sample, sample_errors)
            continue

        # 3. Check if file is empty
        try:
            with gzip.open(fastq_path, "rt") as f:
                first_line = next(f, None)
                if first_line is None:
                    sample_errors.append("FASTQ file is empty.")
                    errors_found = True
                    print_errors(sample, sample_errors)
                    continue
        except Exception as e:
            sample_errors.append(f"Error reading FASTQ file: {e}")
            errors_found = True
            print_errors(sample, sample_errors)
            continue

        # 4. Check if file is in FASTQ format (first few records)
        try:
            with gzip.open(fastq_path, "rt") as f:
                lines = [next(f) for _ in range(8)]
                for i in range(0, len(lines), 4):
                    if not lines[i].startswith("@"):
                        sample_errors.append(
                            "Malformed FASTQ: Read does not start with '@'"
                        )
                        break
        except Exception as e:
            sample_errors.append(f"Error parsing FASTQ format: {e}")

        # 5. Check sample name against wildcard constraint
        if not sample_pattern.fullmatch(sample):
            sample_errors.append(
                "Sample name violates wildcard constraint: check README.md"
            )

        if sample_errors:
            errors_found = True
            print_errors(sample, sample_errors)

    if errors_found:
        print(
            "\nERROR: One or more input validation errors occurred. Aborting Snakemake run.\n"
        )
        sys.exit(1)
    else:
        print("\nInput validation passed!\n")


def print_errors(sample, errors):
    print(f"\nERROR: Sample '{sample}' failed validation with the following issues:")
    for err in errors:
        print(f"   - {err}")
