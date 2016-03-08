# ----------------------------------------------------
# MODULES
# ----------------------------------------------------

from Bio import SeqIO
import argparse
import os
import sys
import socket
import time
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
from Bio.Blast import NCBIXML
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO



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
def run_blast(in_file, db, verbose):
    '''
    This function runs BLAST and prints FASTA files with orthologs
    for each input protein.
    '''

    # GET SEQUENCE DICTIONARIES
    query_dict  = fasta_to_dict(in_file, verbose)
    target_dict = fasta_to_dict(db + ".fa", verbose)

    # RUN BLAST
    if verbose:
        sys.stderr.write("# Running NCBI-BLASTp\n")
    blastp_cmd = Blastp(
        query    = in_file,
        db       = db,
        evalue   = 0.001,
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
        out = open("tmp/" + str(i) + ".fa", "w")
        query_seq = query_dict[blast_record.query]
        out.write(">" + blast_record.query + "\n" + query_seq + "\n")

        for alignment in blast_record.alignments[0:10]:
            target_seq = target_dict[alignment.hit_def]
            best_eval = alignment.hsps[0].expect
            out.write(">" + alignment.hit_def + "\n" + target_seq + "\n")
            for hsp in alignment.hsps:
                if hsp.expect <= best_eval:
                    best_eval = hsp.expect
            print(best_eval)


# ----------------------------------------------------
def fasta_to_dict(fasta, verbose):
    '''
    Creates a dictionary with the FASTA headers as keys and their
    sequences as values.
    '''
    handle      = open(fasta, "rU")
    record_dict = dict()
    if verbose:
        sys.stderr.write("# Reading FASTA file %s\n" % fasta)

    for record in SeqIO.parse(handle, "fasta"):
        record_dict[record.description] = str(record.seq)
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
def do_msa(fasta):
    '''
    Creates a MSA with all the sequences in the fasta file
    '''

# ----------------------------------------------------
# MAIN
# ----------------------------------------------------

if options.verbose:
    print_start_rep()

create_directories()

run_blast(options.input, options.database, options.verbose)



first = open("tmp/1.fa", "r")
do_msa(first)
second = open("tmp/2.fa", "r")
do_msa(second)
