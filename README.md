# MIRROR TREE

## Creating BLAST database

> Run these commands in your working directory

```sh
mkdir db
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -O db/uniprot_sprot.dat.gz
zcat db/uniprot_sprot.dat.gz > db/all_uniprot.fa
makeblastdb -in db/all_uniprot.fa -out db/all_uniprot -dbtype prot
```
