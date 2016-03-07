# ----------------------------------------------------
# MODULES
# ----------------------------------------------------

from Bio import SeqIO
import argparse


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

options = parser.parse_args()


# ----------------------------------------------------
# MAIN
# ----------------------------------------------------

# READ INPUT FASTA FILE
fh = open(options.input, "rU")
sequence_dict = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
fh.close()
