from Bio import SeqIO
import argparse
import sys


parser = argparse.ArgumentParser(description="Reads SWISSPROT and writes FASTA.")
parser.add_argument(
    "-i", "--input",
    dest     = "input",
    action   = "store",
    default  = None,
    required = True,
    help     = "Input SWISSPROT file."
)

parser.add_argument(
    "-o", "--output",
    dest     = "output",
    action   = "store",
    default  = None,
    required = True,
    help     = "Output FASTA file."
)

options = parser.parse_args()

# READ SWISSPROT FILE AND SAVE IT AS FASTA
input_handle  = open(options.input, "rU")
output_handle = open(options.output, "w")

sequences = SeqIO.parse(input_handle, "swiss")
count     = SeqIO.write(sequences, output_handle, "fasta")

output_handle.close()
input_handle.close()

sys.stderr.write("Converted %i records\n" % count)
