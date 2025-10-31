import os
import json
import shutil


def concatenate_logs(log_dir, output):
    logs = {}
    for rule in sorted(os.listdir(log_dir)):
        rule_path = os.path.join(log_dir, rule)
        if not os.path.isdir(rule_path):
            continue

        rule_logs = {}
        for log_file in sorted(os.listdir(rule_path)):
            file_path = os.path.join(rule_path, log_file)

            # Skip empty files
            if os.path.getsize(file_path) == 0:
                continue

            # Use filename without extension as sample name
            sample = os.path.splitext(log_file)[0]

            # Read log text
            with open(file_path, "r", encoding="utf-8", errors="replace") as f:
                text = f.read().strip()

            rule_logs[sample] = text

        if rule_logs:  # only keep non-empty rules
            logs[rule] = rule_logs

    # Write JSON
    with open(os.path.join(output), "w", encoding="utf-8") as out:
        json.dump(logs, out, indent=2, ensure_ascii=False)

    print("logs.json written with", sum(len(v) for v in logs.values()), "log entries.")

    # Remove all subfolders in the log directory
    for rule in os.listdir(log_dir):
        rule_path = os.path.join(log_dir, rule)
        if os.path.isdir(rule_path):
            shutil.rmtree(rule_path)

    print("All subfolders removed, only logs.json remains.")


if __name__ == "__main__":
    log_dir = "logs/"
    output_file = str(snakemake.output)
    concatenate_logs(log_dir, output_file)
