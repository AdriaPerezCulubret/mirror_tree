#!/usr/bin/python3
# ----------------------------------------------------
# MODULES
# ----------------------------------------------------

import os
import re
import sys
import time
import glob
import gzip
import socket
import argparse
import warnings
import itertools
import subprocess
import Mascarpone
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO


# ----------------------------------------------------
# OPTIONS
# ----------------------------------------------------

parser = argparse.ArgumentParser(description="Predicts protein-protein interactions.")
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
    help     = "FASTA file with sequences to find homologs."
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

parser.add_argument(
    "-ints", "--ints",
    dest     = "ints",
    help     = "Tabular file with interactions to consider. Used for train/testing."
)

parser.add_argument(
    "-t", "--taxonomy",
    dest     = "taxonomy",
    action   = "store",
    default  = None,
    required = True,
    help     = "Tabular file with species."
)

parser.add_argument(
    "-c", "--cutoff",
    dest     = "cutoff",
    type     = float,
    default  = 0.5,
    help     = "r cutoff to classify interacting protein pairs."
)

parser.add_argument(
    "-R", "--RNA",
    dest     = "RNA",
    action   = "store",
    default  = None,
    required = True,
    help     = "FASTA file (gzipped) from Silva with aligned rRNA."
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
def print_job(string):
    sys.stderr.write('''
# ------------------------------------------
# %s
# ------------------------------------------
''' % string)


# ----------------------------------------------------
def test_program(program):
    try:
        # pipe output to /dev/null for silence
        null = open("/dev/null", "w")
        subprocess.Popen(program, stdout=null, stderr=null)
        null.close()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            raise Mascarpone.ProgramNotFound(program)
        else:
            raise Exception("Ooops...Something went wrong")

# ----------------------------------------------------
def test_all():
    test_program("jackhmmer")
    test_program("hmmscan")
    test_program("hmmfetch")
    test_program("hmmalign")

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
        sys.stderr.write("# Reading FASTA file %s..." % fasta)
        sys.stderr.flush()

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
        record_dict[record.id] = obj

    handle.close()

    if verbose:
        sys.stderr.write("ok\n")
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
def share_homolog_sp(seq1, seq2, k, tax):
    '''
    Checks if two seq objects share at least k homologs in k species
    '''
    common = seq1.get_homolog_species() & seq2.get_homolog_species() & tax
    if len(common)>= k:
        return common
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

    for seq_hom_obj in sorted(seqobj.homologs, key = lambda seq : seq.species ):
        if seq_hom_obj.species in common_sp:
            SeqIO.write(seq_hom_obj, output_handle, "fasta")
        else:
            continue
    output_handle.close()

    return

# ----------------------------------------------------
def erase_temp(verbose):
    '''
    Removes temporary files
    '''
    if verbose:
        sys.stderr.write("# Removing all tmp files...\n")
    files = glob.glob('tmp/*')
    for f in files:
        os.remove(f)
    if verbose:
        sys.stderr.write("# Removed all tmp files\n\n#  Bye!\n")

# ----------------------------------------------------
def read_interactome(int_file):
    interactions = set()
    intis = open(int_file, "r")
    next(intis) # Skip first line
    exp = re.compile('(uniprotkb:)?(.+)')
    for line in intis:
        lista = line.split("\t")
        try:
            matchA = exp.search(lista[0]).groups()[1]
            matchB = exp.search(lista[1]).groups()[1]
            interactions.add((matchA,matchB))
            interactions.add((matchB,matchA))
        except:
            continue

    return interactions

# ----------------------------------------------------

def read_jack(file, query_dict, target_dict):
    fh = open(file, "r")

    already_added_sp = set()
    for line in fh:
        if line[0] == "#":
            continue
        cols    = line.split()
        t, q, e = cols[0], cols[2], cols[4]
        if (q, target_dict[t].species) not in already_added_sp:
            query_dict[q].homologs.append(target_dict[t])
            already_added_sp.add((q,target_dict[t].species))
        else:
            continue

    return query_dict

# ----------------------------------------------------
def run_jackhmmer(query_dict, target_dict, query, target, verbose):
    '''
    '''
    if verbose:
        sys.stderr.write("# Performing jackhmmer search...")
        sys.stderr.flush()
    #os.system("jackhmmer -E 1e-5 -N 3 --tblout tmp/jackhmmer.tbl --chkhmm tmp/chkhmm %s %s > /dev/null" % (query, target))
    query_dict = read_jack("tmp/jackhmmer.tbl", query_dict, target_dict)
    if verbose:
        sys.stderr.write("ok\n")
        sys.stderr.flush()

    return query_dict

# ----------------------------------------------------
def hmmer_align(seq_file, hmm_file):
    '''
    Aligns all the sequences homologous to the HMM generated by jackhmmer
    '''
    cmd = 'hmmalign --outformat Stockholm "%s" "%s" | perl -ne \'print uc($_);\' > %s.out' %(hmm_file, seq_file, seq_file)
    os.system(cmd)

# ----------------------------------------------------
def hmm_fetch_finder(qname, clean_name):
    '''
    This function searches the HMM for the given sequence name iteratively in the different
    chkhmm files generated by jackhmmer
    '''
    for i in reversed( range(3) ): # Because N = 3 for JACKHMMER
        filenum    = i + 1
        model_name = qname
        filename   = "tmp/chkhmm-%s.hmm" % filenum

        if filenum != 1:
            model_name = qname + "-" + str(i)

        err_code = os.system("hmmfetch %s '%s' > tmp/%s.hmm 2> /dev/null" % (filename, model_name, clean_name))
        if err_code == 0:
            # HMM was found
            break

# ----------------------------------------------------
def print_hmm(query_dict, verbose):
    '''
    Fetches hmm file from jackhmmer and prints the HMM to a file
    '''
    if verbose:
        sys.stderr.write("# Fetching hmm profiles...")
        sys.stderr.flush()

    for qname, qobj in query_dict.items():
        clean_name = qname.replace("|", "_")
        hmm_fetch_finder(qname, clean_name)

    if verbose:
        sys.stderr.write("ok\n")

    return

# ----------------------------------------------------
def read_taxonomy(file, verbose):
    tax_names = set()
    fh = open(file, "r")
    if verbose:
        sys.stderr.write("# Reading taxonomy file...")
        sys.stderr.flush()
    for line in fh:
        line = line.strip()
        tax_names.add(line)
    if verbose:
        sys.stderr.write("ok\n")
    return tax_names

# ----------------------------------------------------
def print_rRNA(common_sp, rnafile):
    '''
    Prints rRNA seqs for the common species to do an MSA.
    '''



# ----------------------------------------------------
# MAIN
# ----------------------------------------------------

# STARTING PROGRAM
if options.verbose:
    print_start_rep()
create_directories()

# TESTING INSTALLED PROGRAMS
test_all()

# READ PROBLEM SEQUENCES
if options.verbose:
    print_job("READING FASTA FILES")
query_dict  = fasta_to_dict(options.input,    options.verbose)
target_dict = fasta_to_dict(options.database, options.verbose)

# READ TAXONOMY FILE
if options.verbose:
    print_job("READING TAXONOMY FILE")
tax_names = read_taxonomy(options.taxonomy, options.verbose)

# HMMER (JACKHMMER + HMMFETCH)
if options.verbose:
    print_job("DOING HMMER SEARCHES")
query_dict  = run_jackhmmer(
    query_dict,
    target_dict,
    options.input,
    options.database,
    options.verbose
)
#print_hmm(query_dict, options.verbose)

# IF TESTING/TRAINING
if options.ints is not None:
    interactions = read_interactome(options.ints)


# ---------------------------------

# PREDICT INTERACTIONS
i = 1
for seq in itertools.combinations(query_dict.keys(), 2):
    seq1, seq2 = query_dict[ seq[0] ], query_dict[seq[1]]
    seqfile_1 = "tmp/%s_1MSA.fa" % i
    seqfile_2 = "tmp/%s_2MSA.fa" % i
    hmmfile_1 = "tmp/%s.hmm" % seq1.id.replace("|", "_")
    hmmfile_2 = "tmp/%s.hmm" % seq2.id.replace("|", "_")

    # IF TRAINING/TESTING
    if options.ints is not None:
        seq1_id = seq1.id.split("|")[1]
        seq2_id = seq2.id.split("|")[1]
        tup = (seq1_id, seq2_id)
        if tup not in interactions:
            continue

    if options.verbose:
        print_job("INTERACTION: %s <-> %s" % (seq1.id, seq2.id) )

    # Check if they share at least K species from tax_names
    common_sp = share_homolog_sp(seq1, seq2, options.species, tax_names)

    if common_sp is None:
        if options.verbose:
            sys.stderr.write("# They don't have %s common species.\n" % options.species)
        continue

    print_seqs_MSA(seq1, seqfile_1, common_sp)
    print_seqs_MSA(seq2, seqfile_2, common_sp)
    print_rRNA(common_sp, options.RNA)

    if options.verbose:
        sys.stderr.write("# Performing MSA for %s\n" % seq1.id )
    hmmer_align(seqfile_1, hmmfile_1)

    if options.verbose:
        sys.stderr.write("# Performing MSA for %s\n" % seq2.id )
    hmmer_align(seqfile_2, hmmfile_2)

    interaction = Mascarpone.Interaction(seq1, seq2)
    interaction.set_dist_matrix(1, seqfile_1 + ".out")
    interaction.set_dist_matrix(2, seqfile_2 + ".out")
    i += 1
    sys.stdout.write(seq1.id + " ")
    sys.stdout.write(seq2.id + " ")
    sys.stdout.write(str(interaction.get_corr()[0]) + " " + str(interaction.get_corr()[1]) + "\n")
    sys.stdout.flush()

print_job("REMOVING TEMP FILES")
#erase_temp(options.verbose)
