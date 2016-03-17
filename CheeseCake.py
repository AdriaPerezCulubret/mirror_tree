#!/usr/bin/python3
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
import warnings
import itertools
import subprocess
import Mascarpone
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO


# ----------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------
def parse_options():

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
        "-o", "--output",
        dest     = "output",
        action   = "store",
        default  = None,
        required = True,
        help     = "Output filename."
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
        default  = 10,
        help     = "Minimum number of common species to create mirror tree."
    )

    parser.add_argument(
        "-ints", "--ints",
        dest     = "ints",
        help     = "Tabular file with interactions to consider. Used for training/testing. Only for debugging."
    )

    parser.add_argument(
        "-d", "--data",
        dest     = "data",
        action   = "store",
        default  = None,
        required = True,
        help     = "Data directory."
    )

    parser.add_argument(
        "-c", "--cutoff",
        dest     = "cutoff",
        type     = float,
        default  = 0.5,
        help     = "r cutoff to classify interacting protein pairs."
    )

    options = parser.parse_args()
    return options

'''
Prints a start log
'''
def print_report():
    HOSTNAME = socket.gethostname()
    sys.stderr.write("""
# ------------------------------------------
#                  CheeseCake
#                  host: %s
#                  time: %s
# ------------------------------------------
""" %(HOSTNAME, time.ctime()))


'''
Prints a Job log
'''
def print_job(string):
    sys.stderr.write('''
# ------------------------------------------
# %s
# ------------------------------------------
''' % string)


'''
Tests if a given program is installed
'''
def test_program(program):
    try:
        # Pipe output to /dev/null for silence
        null = open("/dev/null", "w")
        subprocess.Popen(program, stdout=null, stderr=null)
        null.close()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            raise Mascarpone.ProgramNotFound(program)
        else:
            raise Exception("Ooops...Something went wrong")

'''
Tests if all the necessary programs are installed
'''
def test_all():
    test_program("jackhmmer")
    test_program("hmmscan")
    test_program("hmmfetch")
    test_program("hmmalign")
    test_program("perl")

'''
Gets the species from a sequence header
'''
def get_species(description):
    p = re.compile('OS=(.*?)\s?(\(.+\))?\s?[\w\d]*?\s*=')
    species = None

    try:
        species = p.search(description).groups()[0]
    except:
        species = "Unknown"

    return species

'''
Creates a dictionary with sequence objects. Keys = sequence ID
'''
def fasta_to_dict(fasta, verbose):
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

'''
Creates all the directories used by the program
'''
def create_directories():
    mypath = "tmp"
    if not os.path.isdir(mypath):
        os.makedirs(mypath)

'''
Checks if two seq objects share at least k homologs in k species
'''
def share_homolog_sp(seq1, seq2, k, tax):
    common = seq1.get_homolog_species() & seq2.get_homolog_species() & tax
    if len(common)>= k:
        return common
    else:
        return None

'''
Creates a FASTA file with the homologs of species common to two
sequences
'''
def print_seqs_MSA(seqobj, filename, common_sp):
    output_handle = open(filename, "w")
    SeqIO.write(seqobj, output_handle, "fasta")

    for seq_hom_obj in sorted(seqobj.homologs, key = lambda seq : seq.species ):
        if seq_hom_obj.species in common_sp:
            SeqIO.write(seq_hom_obj, output_handle, "fasta")
        else:
            continue
    output_handle.close()

    return

'''
Removes temporary files
'''
def erase_temp(verbose):
    if verbose:
        sys.stderr.write("# Removing all tmp files...\n")
    files = glob.glob('tmp/*')
    for f in files:
        os.remove(f)
    if verbose:
        sys.stderr.write("# Removed all tmp files\n\n#  Bye!\n")

'''
Reads the interaction file and returns a set with tuples. Only used for Testing purposes.
User: DO NOT TOUCH
'''
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

'''
Reads jackhmmer tabular output and saves homologs in SequenceObj.
We only keep the first result for each species, as it is the one with the lowest e-value.
'''
def read_jack(file, query_dict, target_dict):
    fh = open(file, "r")
    already_added_sp = set()

    for line in fh:
        if line[0] == "#":
            continue
        cols    = line.split()
        t, q, e = cols[0], cols[2], cols[4]
        if (q, target_dict[t].species) not in already_added_sp:
            if query_dict[q].species != target_dict[t].species:
                # Target species is different from query species
                query_dict[q].homologs.append(target_dict[t])
                already_added_sp.add((q,target_dict[t].species))
        else:
            continue
            
    fh.close()
    return query_dict

'''
Runs jackhmmer and reads it using read_jack()
'''
def run_jackhmmer(query_dict, target_dict, query, target, verbose):
    if verbose:
        sys.stderr.write("# Performing jackhmmer search...")
        sys.stderr.flush()
    #os.system("jackhmmer -E 1e-5 -N 3 --tblout tmp/jackhmmer.tbl --chkhmm tmp/chkhmm %s %s > /dev/null" % (query, target))
    query_dict = read_jack("tmp/jackhmmer.tbl", query_dict, target_dict)
    if verbose:
        sys.stderr.write("ok\n")
        sys.stderr.flush()

    return query_dict

'''
Aligns all the sequences homologous to the HMM generated by jackhmmer
'''
def hmmer_align(seq_file, hmm_file):
    cmd = 'hmmalign --outformat Stockholm "%s" "%s" | perl -ne \'print uc($_);\' > %s.out' %(hmm_file, seq_file, seq_file)
    os.system(cmd)

'''
This function searches the HMM for the given sequence name iteratively in the different
chkhmm files generated by jackhmmer
'''
def hmm_fetch_finder(qname, clean_name):
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

'''
Fetches hmm file from jackhmmer and prints the HMM to a file
'''
def print_hmm(query_dict, verbose):
    if verbose:
        sys.stderr.write("# Fetching hmm profiles...")
        sys.stderr.flush()

    for qname, qobj in query_dict.items():
        clean_name = qname.replace("|", "_")
        hmm_fetch_finder(qname, clean_name)

    if verbose:
        sys.stderr.write("ok\n")

    return

'''
Reads the taxonomy files to build the trees using only these species
'''
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

'''
Prints rRNA seq alignments for the common species
'''
def print_rRNA(common_sp, rnafile, filenum, species):
    align = AlignIO.read(open(rnafile), "stockholm")
    common = common_sp
    common.add(species)
    seqs = list()
    for record in align:
        sp = record.id.replace("_", " ")
        if sp in common:
            seqs.append(record)
        else:
            continue

    if len(seqs) == len(common_sp):
        outname = "tmp/%s.sto" % filenum
        outfh = open(outname, "w")
        seqs = sorted(seqs, key = lambda seq : seq.id.replace("_", " "))
        SeqIO.write(seqs, outfh, "stockholm")
        outfh.close()
    else:
        raise Exception("Something wrong: rRNA length != common_sp\n")

'''
Gets the necessary files in the data directory
'''
def get_data_files(direct):
    taxonomy = direct + "/animals.restricted"
    RNAfile  = direct + "/animals_align.sto"
    return taxonomy, RNAfile

'''
Creates an interaction object and sets the necessary distance matrices.
'''
def define_interacton(seq1, seq2, seqfile_1, seqfile_2, i):
    interaction = Mascarpone.Interaction(seq1, seq2)
    interaction.set_dist_matrix(1, seqfile_1 + ".out")
    interaction.set_dist_matrix(2, seqfile_2 + ".out")
    interaction.set_dist_matrix(3, "tmp/" + str(i) + ".sto")
    return interaction

'''
Writes the output to a given file
'''
def write_output(id1, id2, interaction, cutoff, outfh):
    pears, spear, radj  = interaction.get_corr()
    outfh.write("%s\t%s\t%s\t%s\t%s\n" %(id1, id2, pears, spear, radj))
    return

'''
Performs MSA using hmmalign
'''
def hmmalign(sfile1, sfile2, hfile1, hfile2):
    hmmer_align(sfile1, hfile1)
    hmmer_align(sfile2, hfile2)
    return

'''
Computes the correlations for each possible interaction pair and prints them
'''
def predict_interactions(query_dict, options, RNAfile, tax_names, outfh, interactions):
    i = 1
    for seq in itertools.combinations(query_dict.keys(), 2):
        seq1, seq2 = query_dict[ seq[0] ], query_dict[seq[1]]

        if seq1.species != seq2.species:
            continue # Two sequences from different species can't interact

        seqfile_1 = "tmp/%s_1MSA.fa" % i
        seqfile_2 = "tmp/%s_2MSA.fa" % i
        hmmfile_1 = "tmp/%s.hmm" % seq1.id.replace("|", "_")
        hmmfile_2 = "tmp/%s.hmm" % seq2.id.replace("|", "_")

        # IF TRAINING/TESTING
        if options.ints is not None:
            # Assumes that gene_name is like: sp|NAME|something
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

        # PRINT SEQS FOR MSA
        print_seqs_MSA(seq1, seqfile_1, common_sp)
        print_seqs_MSA(seq2, seqfile_2, common_sp)
        print_rRNA(common_sp, RNAfile, i, seq1.species)

        # HMM ALIGN
        if options.verbose:
            sys.stderr.write("# Performing MSA for %s\n" % seq1.id )
        hmmalign(seqfile_1, seqfile_2, hmmfile_1, hmmfile_2)

        interaction = define_interacton(seq1, seq2, seqfile_1, seqfile_2, i)
        write_output(seq1.id, seq2.id, interaction, options.cutoff, outfh)
        i += 1

'''
MAIN PROGRAM
'''
def main():
    # GET THE OPTIONS
    options = parse_options()

    # STARTING PROGRAM
    if options.verbose:
        print_report()
    create_directories()

    # TESTING INSTALLED PROGRAMS
    test_all()
    taxonomy, RNAfile = get_data_files(options.data)

    # READ PROBLEM SEQUENCES
    if options.verbose:
        print_job("READING FASTA FILES")
    query_dict  = fasta_to_dict(options.input,    options.verbose)
    target_dict = fasta_to_dict(options.database, options.verbose)

    # READ TAXONOMY FILE
    if options.verbose:
        print_job("READING TAXONOMY FILE")
    tax_names = read_taxonomy(taxonomy, options.verbose)

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
    interactions = None
    if options.ints is not None:
        interactions = read_interactome(options.ints)

    # OUTPUT
    outfh = open(options.output, "w")
    outfh.write("SEQ1\tSEQ2\tPearson\tSpearman\tr_Adjusted\n")
    predict_interactions(query_dict, options, RNAfile, tax_names, outfh, interactions)
    outfh.close()

    print_job("REMOVING TEMP FILES")
    #erase_temp(options.verbose)
    if options.verbose:
        print_report()


main()
