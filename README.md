# MIRROR TREE

## Creating BLAST database

> Run these commands in your working directory

```sh
mkdir db
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
cd db/
tar -zxvf uniprot_sprot_human.dat.gz
cd ..
python3 bin/db_generator.py -i db/uniprot_sprot_human.dat -o db/human_uniprot.fa
makeblastdb -in db/human_uniprot.fa -out db/human_uniprot -dbtype prot
```
