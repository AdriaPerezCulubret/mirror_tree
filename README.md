# MIRROR TREE

## Creating BLAST database

> Run these commands in your working directory

```sh
mkdir db
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -O db/uniprot_sprot.dat.gz
zcat db/uniprot_sprot.dat.gz > db/all_uniprot.fa
makeblastdb -in db/all_uniprot.fa -out db/all_uniprot -dbtype prot
```


> Get the databases
```sh
wget 'ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt'
python3 create_sets.py -i intact.txt -db db/uniprot_sprot.fasta -o input.fasta
```

create_sets fetch circa than 77.000 IDs but gets circa 45.000 Ids from our uniprot_swissprot.fasta file, be aware. 




```r
hola <- read.table(file="prueba")
library(ggplot2)
ggplot(hola) + geom_point(aes(x=V4, y=V3), position="jitter")
```

```sh
wget 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam29.0/Pfam-A.hmm.gz'
```