# MIRROR TREE

## Creating BLAST database

```sh
mkdir db
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_human.dat.gz -O db/uniprot_sprot_human.dat.gz
cd db/
tar -zxvf uniprot_sprot_human.dat.gz
cd ..
python3 bin/db_generator.py -i db/uniprot_sprot_human.dat -o db/human_uniprot.fa
```
