import os
import tarfile
import shutil
import sys

def tar_and_cleanup(results_dir, output_tar, exclude_dirs):
    # Open tarball for writing
    with tarfile.open(output_tar, "w:gz") as tar:
        for item in os.listdir(results_dir):
            item_path = os.path.join(results_dir, item)
            if item in exclude_dirs:
                print(f"Skipping excluded: {item}", file=sys.stderr)
                continue
            if os.path.islink(item_path):
                print(f"Skipping symlink: {item}", file=sys.stderr)
                continue
            print(f"Adding to tar: {item}", file=sys.stderr)
            tar.add(item_path, arcname=item)

    # Remove all non-excluded items after archiving
    for item in os.listdir(results_dir):
        item_path = os.path.join(results_dir, item)
        if item in exclude_dirs:
            continue
        if os.path.isdir(item_path):
            shutil.rmtree(item_path)
        else:
            os.remove(item_path)

if __name__ == "__main__":
    results_dir = "results"
    output_tar = "results/single_sample_similarity_search_data.tar.gz"
    exclude_dirs = ["abundance", "boxplots", "qc", "test_report.zip","single_sample_similarity_search_data.tar.gz"]

    tar_and_cleanup(results_dir, output_tar, exclude_dirs)
