Attrition Plot

The Plot shows the number for:

- all fasta files (Number of FastQ input reads), 

- the sum of generated hits through similarity search before (Merged similarity hits),

- and after filtering (Filtered fusion reads)

In front of the second plot, there is a thinner stacked bar showing the count discrepancy between unfiltered and filtered hits.
The reasons for filtering are explained as following:
- "Diamond hits < similarity threshold": All hits that have a smaller query/subject - percentage identity as the minimum threshold defined in the config file

- "Diamond hits NOT highest percentage identity per query":  Only the hits with the maximum percentage identity per query ID are accepted. The number of hits filtered for this reason can be seen here

- "Usearch hits > similarity threshold": analogous to "ABR < similarity threshold"

- "Usearch hits NOT highest percentage identity per query": analogous to "ABR Hit not max identity for query ID"

- "Query hit in only one of two databases": In the final filtering step, only those hits gets accepted for which querys can be found in both databases after the previous filterings

Raw data used for this plot can be found in "4. QC/Count Overview Table"