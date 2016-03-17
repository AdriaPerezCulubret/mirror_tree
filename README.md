# MIRROR TREE
			A program created by Sergio Castillo, Joan Martí & Adrià Pérez

## Creating databases

>>> Run these commands in your working directory<<<

> Create uniprot_reposit
```sh
mkdir db
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -O db/uniprot_sprot.dat.gz
zcat db/uniprot_sprot.dat.gz > db/all_uniprot.fa
```

> Get FASTA sequences from proteins that interact from IntAct (BioGrid)
```sh
wget ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt -O db/intact.txt
python3 bin/create_sets.py -i db/intact.txt -db db/uniprot_sprot.fasta -o input.fasta
```

> Get pfam database
```sh
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam29.0/Pfam-A.hmm.gz -O db/Pfam-A.hmm.gz
hmmpress db db/Pfam-A.hmm
hmmscan db/Pfam-A.hmm prova.fa
```

> Get a list of non-interacting proteins as a true negative from Negatome databse
```sh 
wget http://mips.helmholtz-muenchen.de/proj/ppi/negatome/manual_stringent.txt -O db/NonInteract.tbl
python3 bin/create_sets.py -i db/NonInteract.tbl  -db db/uniprot_sprot.fasta -o NonInteract.fasta
```

> Get the alignment of 16S rRNAs for enhancing the correlations
```sh
wget http://www.arb-silva.de/fileadmin/silva_databases/release_123/Exports/SILVA_123_LSURef_tax_silva_full_align_trunc.fasta.gz -O db/SILVA_123_LSURef_tax_silva_full_align_trunc.fasta.gz
```

## EVALUATION



## USAGE
```sh
python3 CheeseCake.py -i FASTA QUERY FILE -db db/all_uniprot.fa -v -sp 10 -t species/animals.tbl -ints LIST INTERACTIONS CHECK OUT > prova15_animalsnoint.out
```
