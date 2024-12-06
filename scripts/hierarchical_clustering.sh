#!/bin/bash

DATA_FOLDER=./data
SCRIPTS_FOLDER=./scripts

mkdir $DATA_FOLDER/clustering

bash $SCRIPTS_FOLDER/dealign.sh $DATA_FOLDER/new_aln.faa > $DATA_FOLDER/clustering/new_dealn.faa

cd-hit -i $DATA_FOLDER/clustering/new_dealn.faa -c 0.9 -o $DATA_FOLDER/clustering/ch_run1.faa
cd-hit -i $DATA_FOLDER/clustering/ch_run1.faa -c 0.7 -o $DATA_FOLDER/clustering/ch_run2.faa
cd-hit -i $DATA_FOLDER/clustering/ch_run2.faa -c 0.5 -n 3 -o $DATA_FOLDER/clustering/ch_run3.faa
cd-hit -i $DATA_FOLDER/clustering/ch_run3.faa -c 0.4 -n 2 -o $DATA_FOLDER/clustering/ch_run4.faa
