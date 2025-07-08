# .snakedeploy/post-deploy.py

import yaml
import os

# Resolve paths relative to the deployment root
base_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(base_dir, "..", "config", "config.yaml")
fastq_dir = os.path.normpath(os.path.join(base_dir, "..", "data", "fastq"))

# 1. Patch config.yaml
try:
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    config["input_validation"] = "no"

    with open(config_path, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)

    print(f"✅ Patched config.yaml: Set input_validation to 'no'.")

except Exception as e:
    print(f"⚠️ Could not patch config.yaml: {e}")

# 2. Ensure fastq directory exists
try:
    os.makedirs(fastq_dir, exist_ok=True)
    print(f"📁 Ensured directory exists: {fastq_dir}")
except Exception as e:
    print(f"⚠️ Could not create data/fastq directory: {e}")
