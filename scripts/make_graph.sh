# make a profile from the last level of clustering
makeblastdb -in data/clustering/ch_run4.faa -dbtype prot -out data/clustering/finalpro

# all-vs-all BLAST using a threshold of 10E-42
blastp -db data/clustering/finalpro -query data/clustering/ch_run4.faa -outfmt 6 -out data/clustering/BLASTe42 -num_threads 4 -evalue 10e-42
