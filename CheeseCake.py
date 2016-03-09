# ----------------------------------------------------
# MODULES
# ----------------------------------------------------

import argparse
import os
import sys
import glob
import pprint
import socket
import itertools
import time
import re
from Bio import SeqIO
import Mascarpone
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
    dest     = "input",
    action   = "store",
    default  = None,
    required = True,
    help     = "Input FASTA file."
)
parser.add_argument(
    "-db", "--database",
    dest     = "database",
    action   = "store",
    default  = None,
    required = True,
    help     = "Database for BLAST."
)
parser.add_argument(
    "-v", "--verbose",
    dest    = "verbose",
    action  = "store_true",
    default = False,
    help    = "Prints log to STDERR."
)
parser.add_argument(
    "-sp", "--species",
    dest     = "species",
    type     = int,
    default  = 5,
    help     = "Minimum number of common species to create mirror tree."
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

    for blast_record in blast_records:
        # We only keep the best hit (e-value) for each target species
        # Because BLAST results are sorted by e-value, we keep the first
        # one we find
        added_sp = set()
        for alignment in blast_record.alignments:
            # Add species homologs
            if target_dict[alignment.hit_def].species not in added_sp:
                query_dict[blast_record.query].homologs[alignment.hit_def] = target_dict[alignment.hit_def]
                added_sp.add(target_dict[alignment.hit_def].species)
            else:
                continue


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

    # Let's create a dictionary using SequenceObj
    for record in SeqIO.parse(handle, "fasta") :
        obj = Mascarpone.SequenceObj(
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
def do_MSA(files, verbose):
    '''
    Creates a MSA with all the sequences in the fasta file
    '''
    for fasta in files:
        if verbose:
            sys.stderr.write("# MSA for %s\n" % fasta)
        tcoffee_cmd = tcoffee(
            infile  = fasta,
            output  = "clustalw",
            outfile = fasta + ".aln",
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
    Checks if two seq objects share at least k homologs in k species
    '''
    if len(seq1.get_homolog_species() & seq2.get_homolog_species()) >= k:
        return seq1.get_homolog_species() & seq2.get_homolog_species()
    else:
        return None

# ----------------------------------------------------
def print_seqs_MSA(seqobj, filename, common_sp):
    '''
    Creates a FASTA file with the homologs of species common to two
    sequences
    '''
    output_handle = open(filename, "w")
    SeqIO.write(seqobj, output_handle, "fasta")

    for seq_hom_name, seq_hom_obj in seqobj.homologs.items():
        if seq_hom_obj.species in common_sp:
            SeqIO.write(seq_hom_obj, output_handle, "fasta")
        else:
            continue
    output_handle.close()

    return


# ----------------------------------------------------
def erase_temp(verbose):
    if verbose:
        sys.stderr.write("# Removing all tmp files...\n")
    files = glob.glob('tmp/*')
    for f in files:
        os.remove(f)
    if verbose:
        sys.stderr.write("# Removed all tmp files\n\n#  Bye!\n")


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
i = 1
for seq in itertools.combinations(query_dict.keys(), 2):
    seq1, seq2 = query_dict[ seq[0] ], query_dict[seq[1]]

    # Now we should run the MSA for each A and B proteins
    # that share at least k species
    common_sp = share_homolog_sp(seq1, seq2, options.species)
    if len(common_sp) >= options.species:
        file_1 = "tmp/%s_1MSA.fa" % i
        file_2 = "tmp/%s_2MSA.fa"  % i
        print_seqs_MSA(seq1, file_1, common_sp)
        print_seqs_MSA(seq2, file_2, common_sp)
        do_MSA((file_1, file_2), options.verbose)
        i += 1
        interaction = Mascarpone.Interaction(seq1, seq2)
        interaction.set_dist_matrix(1, file_1 + ".aln")
        interaction.set_dist_matrix(2, file_2 + ".aln")
        print(str(interaction.get_corr()))





# REMEMBER TO REMOVE ALL THE TMP FILES!!!

erase_temp(options.verbose)

