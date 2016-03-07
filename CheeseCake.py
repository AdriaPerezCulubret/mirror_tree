# ----------------------------------------------------
# MODULES
# ----------------------------------------------------

from Bio import SeqIO
import argparse
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
from Bio.Blast import NCBIXML

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

options = parser.parse_args()


# ----------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------

def run_blast(in_file, db):
    '''
    This function runs BLAST and prints FASTA files with orthologs
    for each input protein.
    '''

    # RUN BLAST
    blastp_cmd = Blastp(
        query    = in_file,
        db       = db,
        evalue   = 0.001,
        outfmt   = 5,
        out      = "tmp/blast_output.xml"
    )

    stdout, stderr = blastp_cmd()

    # PARSE OUTPUT
    result_handle = open("tmp/blast_output.xml")
    blast_records = NCBIXML.parse(result_handle)

    for blast_record in blast_records:
        print(blast_record.query_letters)
        for alignment in blast_record.alignments:
            print ("TITLE: %s\n" % alignment.title)
            for hsp in alignment.hsps:
                print("q: %s t: %s eval: %s\n" % (hsp.query, hsp.sbjct, hsp.expect))


# ----------------------------------------------------
# MAIN
# ----------------------------------------------------

run_blast(options.input, options.database)
