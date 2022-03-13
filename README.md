# IterMegaBLAST
A Sequence Obfuscation Method for Genomic Privacy Protection

# Introduction
With the technological advances in recent decades, determining whole genome sequencing of a person has become feasible and affordable. As a result, large-scale individual genomic sequences are produced and collected for genetic medical diagnoses and cancer drug discovery, which, however, simultaneously poses serious challenges to the protection of personal genomic privacy. It is highly urgent to develop methods which make the personal genomic data both utilizable and confidential. Existing genomic privacy-protection methods are either time-consuming for encryption or with low accuracy of data recovery. To tackle these problems, this paper proposes a sequence similarity-based obfuscation method, namely <b> IterMegaBLAST </b>, for fast and reliable protection of personal genomic privacy.

# Method Description

IterMegaBLAST includes two steps: (1) sequence alignment and clustering; and (2) obfuscation (or anonymization). The whole process for IterMegaBLAST can be summarized as follows: (suppose we used <i>k</i> = 2 for the <i>k</i>-anonymity)
1.	Randomly select a sequence as the initial query sequence from a sequence database with <i>N</i> DNA sequences;
2.	The rest sequences in a dataset is treated as the searching database (<i>N</i>-1);
3.	Use MegaBLAST to generate a list of matched sequences (with high similarity) as well as alignment results between the query and the matched ones;
4.	Select the top similar one as the matched one, and the query and the matched sequences form a cluster (i.e., <i>k</i>=2);
5.	Use the generalization method used in DNALA to produce the generalized sequence of the cluster, and also calculate the distance between the generalized (parent) sequence and the two sequences in the cluster;
6.	Remove these two sequences, treat the rest (<i>N</i>-2) sequences as a new dataset;
7.	Repeat Step 1-6 until all sequences are exhausted;
8.	In case the number of sequences is odd, treat the final three sequences as a cluster. Use MegaBLAST to align these three sequences, and then calculate the average distance among these three distances.

# Usage

We implemented IterMegaBLAST by Perl. To facilitate customizing IterMegaBLAST, we provided all the related Perl scripts. For quick start, users can follow the instructions below:

```{r}
#1. sequence alignment and clustering
./get_size2_cluster.pl.pl <fasta_dir> <tmp dir> <results file>
#2. sequence obfuscation
./get_dist_from_align.pl <ncalign_file> <dist>
```
