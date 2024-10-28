#!/bin/bash

DATA_FOLDER=../data

bash dealign.sh $DATA_FOLDER/new_aln.faa > $DATA_FOLDER/new_dealn.faa

cd-hit -i $DATA_FOLDER/new_dealn.faa -c 0.9 -o $DATA_FOLDER/ch_run1.faa
cd-hit -i $DATA_FOLDER/ch_run1.faa -c 0.7 -o $DATA_FOLDER/ch_run2.faa
cd-hit -i $DATA_FOLDER/ch_run2.faa -c 0.5 -n 3 -o $DATA_FOLDER/ch_run3.faa
cd-hit -i $DATA_FOLDER/ch_run3.faa -c 0.4 -n 2 -o $DATA_FOLDER/ch_run4.faa
