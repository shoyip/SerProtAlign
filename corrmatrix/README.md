These are programs in C++ that allow for the computation of the
$\~{C}_{ij}$ matrix as defined as

$$
\~{C}_{ij} = \sqrt{\sum_{ab} C_{ij}^{ab}^2}
$$

where

$$
C_{ij}^{ab} = f_{ij}^{ab} - f_i^a f_j^b
$$

where again $f_{ij}^{ab}$ is the two-point frequency and
$f_i^a$ is the one-point frequency.

The most efficient implementation by far is the
`corrmatrix_mmap.cpp` version.

Refer to the paper

> Rivoire, Olivier, Kimberly A. Reynolds, and Rama Ranganathan. "Evolution-based functional decomposition of proteins." PLoS computational biology 12.6 (2016): e1004817.
