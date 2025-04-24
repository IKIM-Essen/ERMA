Count Overview

The Plot shows the number for:

- all fasta files (input reads), 

- the sum of generated hits through similarity search before,

- and after filtering

In front of the second plot, there is a thinner stacked bar showing the count discrepancy between unfiltered and filtered hits.
The reasons for filtering are explained as following:
- "ABR < similarity threshold": All hits that have a smaller query/subject - percentage identity as the minimum threshold defined in the config file

- "ABR Hit not max identity for query ID":  Only the hits with the maximum percentage identity per query ID are accepted. The number of hits filtered for this reason can be seen here

- "16S > similarity threshold": analogous to "ABR < similarity threshold"

- "16S Hit not max identity for query ID": analogous to "ABR Hit not max identity for query ID"

- "Query hit not in both databases": In the final filtering step, only those hits gets accepted for which querys can be found in both databases after the previous filterings

Raw data used for this plot can be found in "2. Single Sample Abundance Data/{sample}/Count Overview"