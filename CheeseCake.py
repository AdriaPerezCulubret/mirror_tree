# ----------------------------------------------------
# MODULES
# ----------------------------------------------------

import argparse
import os
import sys
import pprint
import socket
import itertools
import time
import re
from Bio import SeqIO
import BiopythonImproved
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Blast import NCBIXML
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
from Bio.Align.Applications import TCoffeeCommandline     as tcoffee



# ----------------------------------------------------
# OPTIONS
# ----------------------------------------------------

parser = argparse.ArgumentParser(description="OMNOMNOM.")
parser.add_argument(
    "-i", "--input",
    dest="input",
    action="store",
    default=None,
    required=True,
    help="Input FASTA file."
)
parser.add_argument(
    "-db", "--database",
    dest="database",
    action="store",
    default=None,
    required=True,
    help="Database for BLAST."
)
parser.add_argument(
    "-v", "--verbose",
    dest="verbose",
    action="store_true",
    default=False,
    help="Prints log to STDERR."
)


options = parser.parse_args()


# ----------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------

def print_start_rep():
    HOSTNAME = socket.gethostname()
    sys.stderr.write("""
# ------------------------------------------
#                  CheeseCake
#                  host: %s
#            start time: %s
# ------------------------------------------

""" %(HOSTNAME, time.ctime()))

# ----------------------------------------------------
def run_blast(in_file, db, verbose, query_dict):
    '''
    This function runs BLAST and prints FASTA files with orthologs
    for each input protein.
    '''

    # GET SEQUENCE DICTIONARIES
    target_dict = fasta_to_dict(db + ".fa", verbose)

    # RUN BLAST
    if verbose:
        sys.stderr.write("# Running NCBI-BLASTp\n")
    blastp_cmd = Blastp(
        query    = in_file,
        db       = db,
        evalue   = "1e-100",
        outfmt   = 5,
        out      = "tmp/blast_output.xml"
    )

    stdout, stderr = blastp_cmd()

    if verbose:
        sys.stderr.write("# NCBI-BLASTp... ok\n\n")


    # PARSE OUTPUT
    result_handle = open("tmp/blast_output.xml")
    blast_records = NCBIXML.parse(result_handle)

    i = 1
    for blast_record in blast_records:
        out = open("tmp/" + str(i) + ".fa", "a")
        SeqIO.write(query_dict[blast_record.query], out, "fasta")

        for alignment in blast_record.alignments:
            t_species = target_dict[alignment.hit_def].species

            # Add species homologs
            query_dict[blast_record.query].homolog_sp.add(t_species)

            # Add homolog ids with their species
            query_dict[blast_record.query].homologs[t_species] = target_dict[alignment.hit_def].description

            # Write to FASTA
            SeqIO.write(target_dict[alignment.hit_def], out, "fasta")
        i += 1

    return query_dict

# ----------------------------------------------------
def get_species(description):
    '''
    Gets the species from a sequence header
    '''
    p = re.compile('OS=(.*?)\s?(\(.+\))?\s?[\w\d]*?\s*=')
    species = None

    try:
        species = p.search(description).groups()[0]
    except:
        species = None

    return species

# ----------------------------------------------------
def fasta_to_dict(fasta, verbose):
    '''
    Creates a dictionary with sequence objects. Keys = FASTA header
    '''
    handle      = open(fasta, "rU")
    record_dict = dict()
    if verbose:
        sys.stderr.write("# Reading FASTA file %s\n" % fasta)

    record_dict = dict()

    # Let's create a dictionary using SeqRecordOrg
    for record in SeqIO.parse(handle, "fasta") :
        obj = BiopythonImproved.SeqRecordOrg(
            identifier  = str(record.id),
            seq         = record.seq,
            name        = record.name,
            description = record.description,
            dbxrefs     = record.features,
            features    = None,
            annotations = None,
            species     = get_species(record.description),
            letter_annotations = None,
        )
        record_dict[record.description] = obj

    handle.close()

    if verbose:
        sys.stderr.write("# Reading FASTA file %s ... ok\n\n" % fasta)
    return record_dict

# ----------------------------------------------------
def create_directories():
    '''
    Creates all the directories used by the program
    '''
    mypath = "tmp"
    if not os.path.isdir(mypath):
        os.makedirs(mypath)


# ----------------------------------------------------
def do_msa(verbose):
    '''
    Creates a MSA with all the sequences in the fasta file
    '''
    all_files = os.listdir(path="tmp")
    all_files = [ file for file in all_files if file[-2:] == "fa"]


    for fasta in all_files:
        if verbose:
            sys.stderr.write("# MSA for %s\n" % fasta)
        tcoffee_cmd = tcoffee(
            infile  = "tmp/" + fasta,
            output  = "clustalw",
            outfile = "tmp/" + fasta + ".aln",
        )
        stdout, stderr = tcoffee_cmd()
        if verbose:
            sys.stderr.write("# MSA complete for %s ... ok\n\n" % fasta)
    return

# ----------------------------------------------------
def do_filter():
    '''
    Filter the MSA:
        - No paralogous sequences
        - At least 4 sequences from different organisms for each query protein
        - The organisms must be coincident from all proteins
    '''
    all_files = os.listdir(path="tmp")
    all_files = [ file for file in all_files if file[-2:] == "aln"]
    return

# ----------------------------------------------------
def share_homolog_sp(seq1, seq2, k):
    '''
    Checks if two  objects
    '''
    pass


# ----------------------------------------------------
# MAIN
# ----------------------------------------------------

# STARTING PROGRAM
if options.verbose:
    print_start_rep()
create_directories()

# READ PROBLEM SEQUENCES
query_dict  = fasta_to_dict(options.input, options.verbose)

# RUN BLAST
query_dict = run_blast(options.input, options.database, options.verbose, query_dict)

# PREDICT INTERACTIONS
for seq in itertools.combinations(query_dict.keys(), 2):
    seq1, seq2 = query_dict[ seq[0] ], query_dict[seq[1]]

    # Now we should run the MSA for each A and B proteins
    # that share at least k species
    if share_homolog_sp(seq1, seq2, 4):
        print("They share at least 4 species!\n")

# REMEMBER TO REMOVE ALL THE TMP FILES!!!
