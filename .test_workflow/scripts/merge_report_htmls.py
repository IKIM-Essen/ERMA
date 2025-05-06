import glob
import os
import tarfile
import sys


def merge_htmls(path):
    results_dir = path
    os.chdir(results_dir)
    combined_html_path = os.path.join(results_dir, "combined_report.html")
    evalue_boxplot_path = os.path.join(results_dir, "combined_evalue_boxplot.png")
    percidt_boxplot_path = os.path.join(results_dir, "combined_percidt_boxplot.png")
    alignlen_boxplot_path = os.path.join(results_dir, "combined_alignlen_boxplot.png")
    tar_file_path = os.path.join(results_dir, "all_reports.tar.gz")

    sample_files = [file.removesuffix(".fastq.gz") for file in glob.glob("*fastq.gz")]

    with open(combined_html_path, "w") as combined_file:
        combined_file.write(
            "<html><head><title>All Samples Results</title></head><body>"
        )
        combined_file.write("<h1>All Samples Results</h1>")
        combined_file.write(
            "<table style='width:100%; border-collapse: collapse;' border='1'>"
        )

        # Iterate over each sample directory and organize them in a 3x4 matrix
        for i, sample_dir in enumerate(sample_files):
            if i % 1 == 0:
                if i > 0:
                    combined_file.write("</tr>")
                combined_file.write("<tr>")

            sample_path = os.path.join(results_dir, sample_dir)
            if os.path.isdir(sample_path):
                report_path = os.path.join(sample_path, f"{sample_dir}_altair.html")
                if os.path.exists(report_path):
                    combined_file.write("<td style='vertical-align: top; width: 25%;'>")
                    combined_file.write(f"<h2>Results for {sample_dir}</h2>")
                    combined_file.write(
                        f"<iframe src='{os.path.relpath(report_path)}' style='width: 100%; height: 650px; border: none;'></iframe>"
                    )
                    combined_file.write("</td>")

        combined_file.write("</tr>")
        combined_file.write("</table>")
        combined_file.write("</body></html>")

    print("Combined report created at:", combined_html_path)

    with tarfile.open(tar_file_path, "w:gz") as tar:
        # Add the combined report
        tar.add(combined_html_path, arcname=os.path.basename(combined_html_path))

        # Add boxplots
        tar.add(evalue_boxplot_path, arcname=os.path.basename(evalue_boxplot_path))
        tar.add(percidt_boxplot_path, arcname=os.path.basename(percidt_boxplot_path))
        tar.add(percidt_boxplot_path, arcname=os.path.basename(alignlen_boxplot_path))

        # Add each sample's report
        for sample_dir in sample_files:
            sample_path = os.path.join(results_dir, sample_dir)
            if os.path.isdir(sample_path):
                report_path = os.path.join(sample_path, f"{sample_dir}_altair.html")
                if os.path.exists(report_path):
                    tar.add(
                        report_path,
                        arcname=os.path.join(sample_dir, f"{sample_dir}_altair.html"),
                    )

    print("Tar file created at:", tar_file_path)


if __name__ == "__main__":
    path = snakemake.input.path
    sys.stderr = open(snakemake.log[0], "w")
    merge_htmls(path)
