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
from Bio.Blast import NCBIXML
from Bio.SearchIO import HmmerIO
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
        record_dict[record.id] = obj

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
    #SeqIO.write(seqobj, output_handle, "fasta")

    for seq_hom_obj in sorted(seqobj.homologs, key = lambda seq : seq.species ):
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

    for line in fh:
        if line[0] == "#":
            continue
        cols    = line.split()
        t, q, e = cols[0], cols[2], cols[4]
        query_dict[q].homologs.append(target_dict[t])

    return query_dict



# ----------------------------------------------------
def run_jackhmmer(query_dict, target_dict, query, target):
    '''
    '''
    #os.system("jackhmmer -E 1e-20 --tblout tmp/jackhmmer.tbl --chkhmm tmp/chkhmm %s %s > /dev/null" % (query, target))
    query_dict = read_jack("tmp/jackhmmer.tbl", query_dict, target_dict)

    return query_dict

# ----------------------------------------------------
def hmmer_align(seq_file, hmm_file):
    '''
    '''
    os.system("hmmalign --outformat Stockholm -o %s.out %s %s" %(seq_file, hmm_file, seq_file))

# ----------------------------------------------------
def print_hmm(query_dict):
    '''
    Fetches hmm file from jackhmmer and prints the HMM to a file
    '''
    for qname, qobj in query_dict.items():
        clean_name = qname.replace("|", "_")
        os.system("hmmfetch tmp/chkhmm-5.hmm '%s-i4' > tmp/%s.hmm " %(qname, clean_name))

    return



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
query_dict  = fasta_to_dict(options.input,    options.verbose)
target_dict = fasta_to_dict(options.database, options.verbose)

query_dict  = run_jackhmmer(query_dict, target_dict, options.input, options.database)

print_hmm(query_dict)

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
        sys.stderr.write("# Trying to analyze %s and %s...\n" %(seq1.id, seq2.id))
    # Now we should run the MSA for each A and B proteins
    # that share at least k species
    common_sp = share_homolog_sp(seq1, seq2, options.species)
    if common_sp is None:
        if options.verbose:
            sys.stderr.write("# They don't have the necessary common species\n\n")
        continue
    if len(common_sp) >= options.species:
        subset_sp = set()
        for element in range(0, 20):
            try:
                subset_sp.add(common_sp.pop())
            except:
                break
        # Print common species seqs to files
        print_seqs_MSA(seq1, seqfile_1, subset_sp)
        exit(0)
        #print_seqs_MSA(seq2, file_2, subset_sp)

        hmmer_align(seqfile_1, hmmfile_1)

        # MSA!
        #do_MSA((file_1, file_2), options.verbose)

        interaction = Mascarpone.Interaction(seq1, seq2)
        interaction.set_dist_matrix(1, file_1 + ".aln")
        interaction.set_dist_matrix(2, file_2 + ".aln")
        i += 1
        sys.stderr.write(seq1_id + " ")
        sys.stderr.write(seq2_id + " ")
        sys.stderr.write(str(interaction.get_corr()) + "\n")
        sys.stderr.flush()


erase_temp(options.verbose)
