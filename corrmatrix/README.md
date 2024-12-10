These are programs in C++ that allow for the computation of the
$\~{C}_{ij}$ matrix as defined as

$$
\tilde{C}_{ij} = \sqrt{\sum_{ab} C_{ij}^{ab}^2}
$$

where

$$
C_{ij}^{ab} = f_{ij}^{ab} - f_i^a f_j^b
$$

where again $f_{ij}^{ab}$ is the two-point frequency and
$f_i^a$ is the one-point frequency.

The most efficient implementation by far is the
`corrmatrix.cpp` version.

Compile it

```
$ g++ -std=c++20 -fopenmp corrmatrix.cpp -o corrmatrix
```

and use it in the following way

```
$ corrmatrix 369 new_aln.faa weights.csv 10000 corr_matrix.csv
```

where `369` is the number of positions in the alignment, `new_aln.faa` is the
path to the alignment file, `weights.csv` is the path to the file containing
the weights for each sequence, `10000` is the number of subsampled sequences
and `corr_matrix.csv` is the name of the destination file.

Refer to the paper

> Rivoire, Olivier, Kimberly A. Reynolds, and Rama Ranganathan. "Evolution-based functional decomposition of proteins." PLoS computational biology 12.6 (2016): e1004817.
