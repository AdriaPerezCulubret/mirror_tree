from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
import numpy as np

class SequenceObj(SeqIO.SeqRecord):
    '''
    This is a SeqRecord child class with species information
    '''
    def __init__(self, identifier, seq, name, species, description, dbxrefs, features, annotations, letter_annotations):
        super(SequenceObj, self).__init__(
            id = identifier,
            seq = seq,
            name = name,
            description = description,
            dbxrefs = dbxrefs,
            features = features,
            annotations = annotations,
            letter_annotations = letter_annotations
        )
        self.species    = species
        self.homolog_sp = set()
        self.homologs   = dict()

class Interaction(object):
    '''
    This is the interaction Object where it has all the features related with it
    '''

    def __init__(self, seq1, seq2):
        self.matrix1 = list()
        self.matrix2 = list()
        self.seq1 = seq1
        self.seq2 = seq2
        self.correlation = None

    def set_dist_matrix (self, numseq, file):
        aln = AlignIO.read(open(file), 'clustal')
        calculator = DistanceCalculator('blosum62')
        dist_matrix = calculator.get_distance(aln)
        i=0
        j=0
        matrix_list = list()
        for row in dist_matrix:
            j=0
            for column in row:
                if i<j:          # with this, you take out the 0's so n = (N²-N)/2
                    print (dist_matrix[i,j])
                    matrix_list.append(dist_matrix[i,j])
                j+=1
            i+=1
        if numseq == 1:
            self.matrix1 = matrix_list
        elif numseq == 2:
            self.matrix2 = matrix_list
        else:
            raise Exception ("numseq must be either 1 or 2!")

    def get_corr (self):
        if not self.matrix1 or not self.matrix2:
            raise Exception ("Distance matrix not calculated")
        elif self.correlation is None:
            self.correlation = np.corrcoef(self.matrix1,self.matrix2)[1:0]
            return self.correlation
        else: 
            return self.correlation
