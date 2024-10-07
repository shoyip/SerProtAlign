# Serine Protease Alignment

This repository contains the codebase (libraries, codes, notebooks, scripts) for
performing alignments of the Serine Protease family.

## How is the family defined

We can use multiple definitions, they are slightly different between each other:

- The MEROPS database has the S1A Subfamily
- The Pfam database defines the PF00089 Trypsin Family
- We use a manual selection given by Halabi

## The new HMMER alignment

In order to make the new alignment we are going to use a set of bioinformatic
tools.

### The protein database

There are multiple databases of amino acidic sequences. Usually they are created
by feeding all the reads coming from sequencing of genomes, so they are quite
big in size. The most general databases that contain full length sequences of
proteins are:

| Database | Features |
| --- | --- |
| [NCBI Protein Database](https://www.ncbi.nlm.nih.gov/protein) | Contains over 355 milion protein sequence entries, and contains both curated (RefSeq) and non curated sequences. |
| [UniProtKB](https://www.uniprot.org/uniprotkb/) | Contains over 291 milion protein sequences entries, and contains both curated (SwissProt) and non curated (TrEMBL) sequences. |

We are going to choose the UniProt database. It can be downloaded from the
[UniProtKB website](https://www.uniprot.org/uniprotkb), in two separate files
containing SwissProt sequences and TrEMBL sequences. We have just concatenated
the two files in order to get the full UniProtKB database.

> [!WARNING]
> The UniProtKB database amounts to a FASTA file and is huge (the April 2024
> version weighs 112.7 GB). It is not included in this repository but can be
> downloaded by issuing a script.

### The representative MSA

In order to have a good multiple sequence alignment (MSA) we should carefully
select a handful of sequences that we deem to be representative of the family
that we are looking for, and align them. This alignment can be curated either
manually or performed using global alignment algorithms (i.e. MAFFT).

We get a **seed MSA**, which is a small alignment that we think to be good
enough to capture the properties of the target family. We can then proceed to
building an **HMM profile**: HMMs (Hidden Markov Models) are consensus models,
and assume independence of residues between positions. An HMM profile is defined
by the set of transition probabilities from an amino acid state to an insert
state or a delete state, and all the possible combinations among them.

We choose to define our HMM profile from a selection of sequences defined in the
Halabi et al., 2009 paper. As described in the *Experimental Procedures* and in
the *SI Supplemental Experimental Procedures* sections:

- The sequences were obtained by an iterative procedure using PSI-BLAST. This
    means that an initial sequence (say, for example, the UniProtKB `TRY2_RAT`)
    is PSI-BLASTed, yielding to a certain amount of sequence hits. Then they are
    PSI-BLASTed again, and so on and so forth until a convergence is reached
    (i.e. the difference between the sequence set in one run and the next one
    becomes negligible).
- The resulting sequences are full-length and unaligned. They have been aligned
    manually using Cn3D and automatically using ClustalX.
- The alignment has been then truncated to the amino acids present in the
    `TRY2_RAT` sequence, meaning that all positions containing gaps in that very
    sequence were entirely eliminated from the alignment.

In the end this yields to an alignment of 1470 sequences and 223 positions. We
actually use the result of the second step as a seed MSA, which corresponds
hence to the untruncated alignment with 1470 sequences and 832 positions.

We hence run the following command to get the HMM profile, using the `symfrac 0`
option to avoid throwing away gappy positions. In fact, by default `hmmbuild`
sets to 0.5 meaning that it assignes as consensus only columns that include at
most 50% of gaps, while in our case we want to retain every alignment column.

```
$ hmmbuild --symfrac 0 data/serprotalign.hmm data/ref_seqs.faa > logs/hmmbuild_log.out
```

### Search of profile against sequence database

Now that we have a profile we can search for sequences against a sequence
database. This is done by `hmmsearch`: it both looks for sequences that have a
low E-value (i.e. that well match a given HMM profile) and returns the
alignment. *Two birds with one stone!*

```
$ hmmsearch -A data/serprot_matches.stk data/uniprot.faa > logs/hmmsearch_log.out
```

The resulting sequences will be saved to a Stockholm file, which contains many
informations that we do not need in the alignment. So we will just extract the
descriptions and the sequences in a FASTA file. We use [this Perl script](https://biocs-blog.blogspot.com/2010/08/convert-stockholm-sequence-format-to.html)
from Michael Kuhn to perform this operation and to get rid of insertions (dot
and lowercase characters).

```
$ perl scripts/stk2fasta.pl data/serprot_matches.stk > data/serprot_matches.faa
```

### Deduplicating

A drawback of using `hmmsearch` is that it might match the same sequence twice
because two different parts of it had a high enough match with it. Often they
are spurious matches with long sequences, so we get rid of them.

### Throwing away unusual sequences

Before doing any truncation, we would like to find a scheme to reduce the
alignment length without losing too much information. We do this by finding away
sequences that create gaps, and by looking for a threshold such that we don't
loose too many sequences while we gain in compactness of the alignment.

The full alignment has ~10^5 sequences, so we perform this analysis on a
subsample of ~10^4 sequences and apply the resulting statistics to the full
dataset. The analysis is provided in
[a Python Notebook](notebooks/compact_alignment.ipynb).

## The new structural alignment

Another way to perform an alignment is to keep into account the structural
informations. This goes a step further: remember, in fact, that we said that
positions were independent in our previous alignment technique.

We take advantage of the tools provided with the FoldSeek suite in order to
perform this **structural alignment**.

### Choose the seed structures

We choose a handful of structures that we want to be our structural seed. In our
case we choose

| PDB ID | Corresponding UniProt ID | Description |
| --- | --- | --- |
| 1AZZ | P00771 COGS_LEPPG | Brachyurin |
| 1CGH | P08311 CATG_HUMAN | Cathepsin G |
| 1EQ9 | Q7SIG2 CTR1_SOLIN | Chymotrypsin-1 |
| 1EUF | P80219 DDN1_BOVIN | Duodenase-1 |
| 1GDU | P35049 TRYP_FUSOX | Trypsin |
| 1OKX | P00772 CELA1_PIG | Chymotrypsin-like elastase family member 1 |
| 1OSS | P00775 TRYP_STRGR | Trypsin |
| 1T8O | P00766 CTRA_BOVIN | Chymotrypsinogen A |
| 3TGI | P00763 TRY2_RAT | Anionic Trypsin-2 |
| 5DJ7 | Q54137 Q54137_SACER | Trypsin-like protease |

We download their structures in a folder.

### Get the 3di descriptor

The mechanism by which FoldSeek works is that it transforms the amino acid
residues in another type of residues that contain informations about structure
by means of a dictionary (the 3di dictionary). We perform this operation by
issuing the command

```
$ foldseek structureto3didescriptor data/ref_structs/*.pdb data/3di.dump
```

The resulting `3di.dump` file contains many informations: first of all the PDB
structure ID, then the chain identifier, the name of the structure, the amino
acid sequence, the 3di sequence and positional informations about the residues.
We just need the amino acid and 3di sequences for the desired chains.


```
$ awk -F'\t' '{print ">"$1"\n"$3; next}{print}' data/3di.dump > data/3di.faa
```

```
$ awk -F'\t' '{print ">"$1"\n"$2; next}{print}' data/3di.dump > data/3di_aa.faa
```

This means that we have to hand pick the chains of interest, and put the 3di
sequences in one file and the amino acid sequences in another file. We then
align the 3di sequences, then feed the gaps that we get from the 3di alignment
in the amino acid sequences.

```
$ mafft 3di.fasta | awk 'BEGIN {ORS=""} !/^>/ { print $1 } /^>/ {print "\n" $1 "\n" }' | tee data/3di_aligned.faa
```

Now we realign the amino acidic sequences corresponding to the structures by
inserting gaps in the same places as in the 3di sequences.

```
python scripts/realign_structures_aa.py > data/3di_aligned_aa.faa
```

We then use this new alignment as a seed alignment and proceed as described
previously.

## The new SSN-based alignment

TODO

# Prerequisites

Python, Perl, MAFFT, HMMER suite, FoldSeek, wget.

# Impressum

Shoichi Yip. 2024.
