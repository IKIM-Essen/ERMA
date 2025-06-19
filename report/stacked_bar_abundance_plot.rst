Stacked bar abundance plot

The Interactive Boxplot visualizes the distribution of similarity search alignment parameters—alignment length, percentage identity, and e-value—between 16S and antimicrobial resistance (AMR) gene regions across multiple samples.

Each panel in the plot corresponds to one alignment parameter. Within each panel, boxplots compare values for the ABR (AMR gene) and 16S rRNA gene portions across all input samples. This allows users to assess differences in search quality or filtering stringency between gene parts and across conditions or datasets.

Each data point represents a single similarity search hit. Outliers are shown as individual points, and box widths reflect the spread of the data. Axes are scaled linearly for alignment length and percentage identity, and logarithmically for e-values to reflect magnitude differences.

The plot is generated interactively, allowing legend toggling to show or hide gene parts per parameter. Only AMR gene families with at least 1% of the total summed hit counts are included for interpretability.