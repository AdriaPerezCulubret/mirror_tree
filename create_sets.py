# ----------------------------------------------------
# MODULES
# ----------------------------------------------------

import argparse
import os
import sys
import re
from Bio import SeqIO


# ----------------------------------------------------
# OPTIONS
# ----------------------------------------------------

parser = argparse.ArgumentParser(description="OMNOMNOM.")
parser.add_argument(
    "-i", "--interactome",
    dest="inter",
    action="store",
    default=None,
    required=True,
    help="Input interactome."
)
parser.add_argument(
    "-db", "--database",
    dest="database",
    action="store",
    default=None,
    required=True,
    help="Database to retrive from in FASTA format."
)
parser.add_argument(
    "-o", "--output",
    dest="output",
    action="store",
    default=None,
    required=True,
    help="Output file to store the filtered FASTAs."
)

options = parser.parse_args()


# ----------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------

def id_interactome (inter_fh):
    myset = set()
    for line in inter_fh:
        line = line.strip()
        columns = line.split("\t")
        exp = re.compile('uniprotkb:(.+)')
        seqA = columns[0]
        seqB = columns[1]
        matchA = None
        matchB = None
        try:
            matchA = exp.search(seqA).groups()[0]
            matchB = exp.search(seqB).groups()[0]
            myset.add(matchA)
            myset.add(matchB)
        except:
            matchA = None
            matchB = None
    return (myset)


# ----------------------------------------------------

def filter_db (db_fh):
    sys.stdout = open(options.output, 'w')
    db_fh = open (options.database, "rU")
    for entry in SeqIO.parse(db_fh,"fasta"):
        exp = re.compile('sp\|(.*)\|')
        match = None
        try:
            match = exp.search(entry.id).groups()[0]
            if match in prot_inter:
                print (">"+entry.id)
                print (entry.seq)
        except:
            match = None


# ----------------------------------------------------
# MAIN
# ----------------------------------------------------


# Load uniprot ID from verified interactom into a set
inter_fh = open (options.inter)
inter_fh.readline()
prot_inter = id_interactome(inter_fh)
print (len(prot_inter))

# Filter by items on the set in the uniprot data set and print them in a file
filter_db(options.database)




 

