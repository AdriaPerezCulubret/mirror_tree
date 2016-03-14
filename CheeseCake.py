# ----------------------------------------------------
# MODULES
# ----------------------------------------------------

import os
import re
import sys
import time
import glob
import socket
import argparse
import itertools
import subprocess
import Mascarpone
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.SearchIO import HmmerIO

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

parser.add_argument(
    "-ints", "--ints",
    dest     = "ints",
    help     = "Tabular file with interactions to consider. Used for train/testing."
)

parser.add_argument(
    "-p", "--pfam",
    dest     = "pfam",
    action   = "store",
    default  = None,
    required = True,
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
    for line in intis:
        lista = line.split("\t")
        lista[0] = lista[0].split(":")[1]
        lista[1] = lista[1].split(":")[1]
        interactions.add((lista[0],lista[1]))
        interactions.add((lista[1],lista[0]))
    return interactions

# ----------------------------------------------------
def run_hmmscan(query_dict, pfam, input, verbose):
    '''
    This functions searches PFAM domains for our query sequences.
    '''
    os.system("hmmscan -E 1e-20 --domE 1e-10 %s %s > tmp/hmmscan.out" % (pfam, input))
    fh = open("tmp/hmmscan.out", "r")
    hmmer_results = HmmerIO.Hmmer3TextParser(fh)

    for res in hmmer_results:
        query_obj = query_dict[res.id]
        for it in res.hits:
            for hsp in it.hsps:
                dom = Mascarpone.Domain(
                    name   = it.id,
                    evalue = hsp.evalue,
                    coords = (hsp.env_start, hsp.env_end)
                )
                query_obj.domains.append(dom)

    return query_dict

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
    os.system("jackhmmer -E 1e-5 --tblout tmp/jackhmmer.tbl --chkhmm tmp/chkhmm %s %s > /dev/null" % (query, target))
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
    os.system("hmmalign --outformat Stockholm -o %s.out %s %s" %(seq_file, hmm_file, seq_file))
    fh = open("%s.out" % seq_file, "r")
    new = open("KKKK", "w")
    AlignIO.convert(fh, "stockholm", new, "clustal")


# ----------------------------------------------------
def hmm_fetch_finder(qname, clean_name):
    '''
    This function searches the HMM for the given sequence name iteratively in the different
    chkhmm files generated by jackhmmer
    '''
    for i in reversed( range(5) ):
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
def read_mammals(file, verbose):
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
# MAIN
# ----------------------------------------------------

# STARTING PROGRAM
if options.verbose:
    print_start_rep()
create_directories()

# TESTTING INSTALLED PROGRAMS
test_all()

# READ PROBLEM SEQUENCES
if options.verbose:
    print_job("READING FASTA FILES")
query_dict  = fasta_to_dict(options.input,    options.verbose)
target_dict = fasta_to_dict(options.database, options.verbose)

# READ TAXONOMY FILE
if options.verbose:
    print_job("READING TAXONOMY FILE")
tax_names = read_mammals(options.taxonomy, options.verbose)

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
print_hmm(query_dict, options.verbose)

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

    if len(common_sp) >= options.species:
        print_seqs_MSA(seq1, seqfile_1, common_sp)
        print_seqs_MSA(seq2, seqfile_2, common_sp)

        if options.verbose:
            sys.stderr.write("# Performing MSA for %s\n" % seq1.id )
        hmmer_align(seqfile_1, hmmfile_1)

        if options.verbose:
            sys.stderr.write("# Performing MSA for %s\n" % seq2.id )
        hmmer_align(seqfile_2, hmmfile_2)




        #interaction = Mascarpone.Interaction(seq1, seq2)
        #interaction.set_dist_matrix(1, file_1 + ".aln")
        #interaction.set_dist_matrix(2, file_2 + ".aln")
        #i += 1
        #sys.stderr.write(seq1_id + " ")
        #sys.stderr.write(seq2_id + " ")
        #sys.stderr.write(str(interaction.get_corr()) + "\n")
        #sys.stderr.flush()

print_job("REMOVING TEMP FILES")
#erase_temp(options.verbose)
