# MTREE
	A program to predict interactions using mirror trees.
            Created by Sergio Castillo, Joan Martí & Adrià Pérez

## Install mtree
> Move to the root directory of the program and run these commands:

```sh
sudo python3 setup.py install
```


## Usage of mtree

```sh
mtree -i input.fa -db uniprot.fa --data datafolder -c 0.8 -sp 8
```

### Options

* **-h,&nbsp;&nbsp;&nbsp;&nbsp;--help**  

    ```
    Shows the options of the program.
    ```

* **-i,&nbsp;&nbsp;&nbsp;&nbsp;--input** &nbsp;&nbsp;&nbsp;&nbsp; FILENAME

    ```
    REQUIRED. Input FASTA with the sequences to predict interactions.
    ```

* **-o,&nbsp;&nbsp;&nbsp;&nbsp;--output** &nbsp;&nbsp;&nbsp;&nbsp; FILENAME

    ```
    OPTIONAL. Filename of the output. If it's not defined mtree will print to STDOUT.
    ```

* **-db,   --database** &nbsp;&nbsp;&nbsp;&nbsp; FILENAME

    ```
    REQUIRED. Sequence DATABASE in FASTA format. Will be used to find homologs using jackhmmer.
    ```

* **-d,&nbsp;&nbsp;&nbsp;&nbsp;--data** &nbsp;&nbsp;&nbsp;&nbsp; DIRECTORY

    ```
    REQUIRED. Directory with the data files used by mtree.
    ```

* **-sp,   --species** &nbsp;&nbsp;&nbsp;&nbsp; INT

    ```
    OPTIONAL. Number of species needed to build the phylogenetic trees. Default = 10.
    ```

* **-c,&nbsp;&nbsp;&nbsp;&nbsp;--cutoff** &nbsp;&nbsp;&nbsp;&nbsp; FLOAT

    ```
    OPTIONAL. Pearson correlation cutoff to predict interactions. Default = 0.7.
    ```

* **-g,&nbsp;&nbsp;&nbsp;&nbsp;--graph**&nbsp;&nbsp;&nbsp;&nbsp; FILENAME

    ```
    OPTIONAL. If this option is used, mtree will create an html file with the graph drawn using cytoscape.js.
    ```

* **-ints, --ints**&nbsp;&nbsp;&nbsp;&nbsp; FILENAME

    ```
    OPTIONAL. USED ONLY FOR TESTING/DEBUGGING! Considers only interactions in FILENAME.
    ```



### Output of the program
> The program outputs a tabular file to the specified path (using the option -o). This file is of the following form:

|  SEQ1   |  SEQ2   | Pearson   | Spearman   | r_Adjusted   |   Type    |
| :----:  | :----:  | :------:  | :-------:  | :----------: | :-------: |
|  CLOCK  |  BMAL   |   0.7     |   0.6      |     0.8      |  Int      |
|  EGFR   |  TNFA   |   0.3     |   0.2      |     0.4      |  NonInt   |


# Testing the program

## Downloading the databases
> **mtree** needs a database of FASTA sequences to work. We used Swiss-Prot.
> Also, for training purposes, we used INTACT and NEGATOME. These two
> are NOT necessary for the program to work, but we used them to evaluate our method.


### Create directory

```sh
mkdir db
```

### Download uniprot

```sh
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz \
     -O db/uniprot_sprot.dat.gz
```

### Subsetting the database
> Because we are always working with animals we filtered the database. This is NOT required.

#### Downloading the animals list
> Go to the website and download all the animal species names.
> Save it as animals.txt
> This step is NOT necessary. It was only useful for testing.

* http://www.ncbi.nlm.nih.gov/genome/browse/

#### Filter the FASTA
```sh
zcat db/uniprot_sprot.dat.gz | \
perl -e '
%specs = ();
open $fh, "<animals.txt ";
while(<$fh>) {
    chomp;
    $specs{$_} = 1;
};
local $/ = ">";
while(<>){
    chomp;
    ($name, @seq) = split /\n/;
    $sp = "";
    if ($name =~ m/OS=(.*?) GN/) {
        $sp = $1;
    } else {
        next
    };
    if (not exists $specs{$sp}){
        next;
    }  
    print ">$name\n", join("\n", @seq), "\n";
} ' > db/animals_uniprot.fa
```

### Create testing sets
> We needed a set of interacting proteins and a set of non-interacting proteins.
> We downloaded the former from INTACT, and the latter from Negatome

#### Download INTACT
```sh
wget ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt -O db/intact.tbl
```

#### Download NEGATOME

```sh
wget http://mips.helmholtz-muenchen.de/proj/ppi/negatome/manual_stringent.txt -O db/negatome.tbl

# Get some FASTAs from negatome
python3 bin/create_sets.py -i db/NonInt.tbl -db db/animals_uniprot.fa -o 300_negatome.fa --num 300
```





> Get the alignment of 16S rRNAs for enhancing the correlations

```sh
wget http://www.arb-silva.de/fileadmin/silva_databases/release_123/Exports/SILVA_123_LSURef_tax_silva_full_align_trunc.fasta.gz -O db/SILVA_123_LSURef_tax_silva_full_align_trunc.fasta.gz
```




## Plots

```r
library(ggplot2)
data <-read.table(file="testing.tbl")
ggplot(data) +
    geom_point(aes(x=Type, y=Pearson, color=Type), position="jitter") +
    xlab("") + ylab("Coevolution (Pearson)\n") +
    theme_bw() +
    geom_hline(yintercept=0.7, linetype="dashed", alpha=0.6) +
    geom_hline(yintercept=0.8, linetype="dashed", alpha=0.6) +
    geom_hline(yintercept=0.9, linetype="dashed", alpha=0.6) +
    annotate("text", x=0.7, y=0.72, label="P = 0.78") +
    annotate("text", x=0.7, y=0.82, label="P = 0.72") +
    annotate("text", x=0.7, y=0.92, label="P = 0.63") +
    annotate("text", x=2.4, y=0.72, label="R = 0.12") +
    annotate("text", x=2.4, y=0.82, label="R = 0.07") +
    annotate("text", x=2.4, y=0.92, label="R = 0.03")

ggplot(data) +
    geom_point(aes(x=Type, y=Spearman, color=Type), position="jitter") +
    xlab("") + ylab("Coevolution (Spearman)\n") +
    theme_bw() +
    geom_hline(yintercept=0.7, linetype="dashed", alpha=0.6) +
    geom_hline(yintercept=0.8, linetype="dashed", alpha=0.6) +
    geom_hline(yintercept=0.9, linetype="dashed", alpha=0.6) +
    annotate("text", x=0.7, y=0.72, label="P = 0,77") +
    annotate("text", x=0.7, y=0.82, label="P = 0,80") +
    annotate("text", x=0.7, y=0.92, label="P = 0.71") +
    annotate("text", x=2.4, y=0.72, label="R = 0.08") +
    annotate("text", x=2.4, y=0.82, label="R = 0.06") +
    annotate("text", x=2.4, y=0.92, label="R = 0.02")

ggplot(data) +
    geom_point(aes(x=Type, y=r_Adjusted, color=Type), position="jitter") +
    xlab("") + ylab("Coevolution (Partial correlation)\n") +
    theme_bw() +
    geom_hline(yintercept=0.97, linetype="dashed", alpha=0.6) +
    geom_hline(yintercept=0.98, linetype="dashed", alpha=0.6) +
    geom_hline(yintercept=0.99, linetype="dashed", alpha=0.6) +
    annotate("text", x=0.7, y=0.972, label="P = 0.70") +
    annotate("text", x=0.7, y=0.982, label="P = 0.71") +
    annotate("text", x=0.7, y=0.992, label="P = 0.69") +
    annotate("text", x=2.4, y=0.972, label="R = 0.72") +
    annotate("text", x=2.4, y=0.982, label="R = 0.58") +
    annotate("text", x=2.4, y=0.992, label="R = 0.39")
```
