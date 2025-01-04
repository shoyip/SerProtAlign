# Serine Protease Alignment

This repository contains the very first draft of code that was used to prepare a Multiple Sequence Alignment (MSA) of the Serine Protease family.

Data is available on [Zenodo](https://doi.org/10.5281/zenodo.14599586).

## Building an HMM profile

We run the following command to get the HMM profile, using the symfrac 0 option to avoid throwing away gappy positions. In fact, by default hmmbuild sets to 0.5 meaning that it assignes as consensus only columns that include at most 50% of gaps, while in our case we want to retain every alignment column.

```
$ hmmbuild --symfrac 0 serprot.hmm serprot_ref.faa > logs/hmmbuild_log.out
```

## Profile search against a sequence database

We run an HMM search, meaning that the previously obtained HMM profile is used to look up for sequences that are homologous to those in the `serprot_ref.faa` file. This command returns a Stockholm file.

```
$ hmmsearch -A serprot_matches.stk uniprot.faa > logs/hmmsearch_log.out
```

## Converting the MSA from Stockholm format to FASTA and remove inserts

We want a more manageable and slimmer file for the sequence alignment. We get rid of the unwanted parts and of inserts (`.` signs in the alignment) by running [a script written by Michael Kuhn](https://biocs-blog.blogspot.com/2010/08/convert-stockholm-sequence-format-to.html).

```
$ perl scripts/stk2fasta.pl serprot_matches.stk > new_aln0.faa
```

## Cleaning and truncation of the alignment

The alignment is then cleaned and truncated for the purpose of having a better dataset to train Machine Learning models on.
The code is written in Python and is detailed in two Jupyter Notebooks, along with analyses:

- The `notebooks/00_alignment_cleaning.ipynb` notebook contains the code used to go from `new_aln0.faa` to `new_aln.faa`. Sequences containing ambiguity symbols are filtered out, sequences that lack the IVGG pattern are filtered out, columns corresponding to the signal peptide and activation peptides are filtered out, sequences that are length outliers are filtered out.
- The `notebooks/01_compact_alignment.ipynb` notebook contains the code used to go from `new_aln.faa` to `iter_aln.faa`. Three truncation schemes are compared: the **gappy out**, the **compact** and the **compact iterative** scheme. The gappy out scheme is commonly used and it implies the truncation of columns above a certain gap threshold, however it might truncate columns which are relevant for protein function. Instead, we propose a scheme that gets rid of the most amount of columns by leaving the most amount of sequences untrimmed. The iterative scheme improves the statistics of the compact scheme by iteratively applying it for a certain amount of iterations (four in our case).
- The `notebooks/02_alignment_comparison.ipynb` notebook contains the analyses to compare the alignment in `new_aln.faa` and in `iter_aln.faa`. It compares the PCAs of the two alignments, and the analysis of the elements of the correlation matrices of the two alignments.

## Correlation matrix analysis

Correlation matrices contain important informations used by Machine Learning models to learn from biological sequence data.
This is why we compared the distribution of elements of the correlation matrix. To compute the correlation matrix a program in C++ was written, which computes the statistics on a subsample of the alignment. It makes use of the weights computed in the `notebooks/02a_reweighting.ipynb`, that should be run on a machine with GPU.

The program to compute the correlation matrix can be found in the `corrmatrix` folder. It can be run as follows.

```
$ ./corrmatrix 693 new_aln.faa weights.csv 10000 corr_aln.csv
```

In order the parameters are:

1. The length of the alignment (L=693)
2. The file of the alignment (`new_aln.faa`)
3. The file of the weights (`weights.csv`)
4. The number of sequences to be subsampled (M*=10000)
5. The output file (`corr_aln.csv`)

# Impressum

Shoichi Yip. 2025.
