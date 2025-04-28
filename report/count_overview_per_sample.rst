Count Overview per sample

The Table shows an count overview of input, intermediate and output files:

- fasta input: number of reads of the raw fastq/fasta file,

- diamond output: number of hits, the diamond (ABR) similarity search generated on this sample,

- usearch output: number of hits, the usearch (16S) similarity search generated on this sample,

- integration output: number of lines resulting from the merging of diamond and usearch results

- filtered min similarity ABR: lines filtered through user defined minimum percentage identity threshold

- filtered max identity ABR: lines filtered since only the maximum percentage identity per query ID is accepted

- filtered min similarity 16S: lines filtered through user defined minimum percentage identity threshold

- filtered max identity 16S: lines filtered since only the maximum percentage identity per query ID is accepted

- filtered query id mismatch: hits filtered because query ID is not found in both databases

- filtration output: number of lines in the final fight