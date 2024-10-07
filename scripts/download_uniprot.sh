#!/bin/bash

$data = $(dirname "$0")/../data

# download the swissprot and trembl files
wget https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz
wget https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_trembl.fasta.gz

# merge the two files
seqtk mergefa uniprot_sprot.fasta.gz uniprot_trembl.fasta.gz > data/uniprot.faa

# clean the directory
rm uniprot_sprot.fasta.gz
rm uniprot_trembl.fasta.gz
