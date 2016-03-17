import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
import scipy.stats as st
from Bio.Phylo.TreeConstruction import DistanceCalculator

class SequenceObj(SeqIO.SeqRecord):
    '''
    This is a SeqRecord child class with species information.
        species    : The species of the sequence
        homologs   : Dictionary with SequenceObj objects
        homolog_sp : Set with the species of all the homologs of this Sequence
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
        self.species      = species
        self.homologs     = list()
        self.hmm          = None
        self.__homolog_sp = None

    def get_homolog_species(self):
        '''
        This returns a set with the species of all the homologs of this progein
        '''
        if self.__homolog_sp is None:
            self.__homolog_sp = set()
            for homo_seq in self.homologs:
                self.__homolog_sp.add(homo_seq.species)
            return self.__homolog_sp
        else:
            return self.__homolog_sp


class Interaction(object):
    '''
    This is the interaction Object. It has all the attributes/features related with it.
    '''

    def __init__(self, seq1, seq2):
        self.matrix1 = list()
        self.matrix2 = list()
        self.seq1 = seq1
        self.seq2 = seq2
        self.sp_matrix   = None
        self.correlation = None
        self.spearman    = None

    def __matrix_to_list(self,dist_matrix):
        '''
        Transforms a distance matrix to a list.
        '''
        matrix_list = list()
        i=0
        j=0
        for row in dist_matrix:
            j=0
            for column in row:
                if i<j:          # with this, you take out the 0's so n = (N-N)/2
                    matrix_list.append(dist_matrix[i,j])
                j+=1
            i+=1
        return matrix_list

    def __linear_reg_residues(self, x, y):
        '''
        Linear regression that returns the residuals
        '''
        residuals = list()
        A    = np.vstack([x, np.ones(len(x))]).T
        slope, intercept = np.linalg.lstsq(A, y)[0]
        for idx, num in enumerate(x):
            predict = slope * num + intercept
            res     = y[idx] - predict
            residuals.append(res)
        residuals = np.append(residuals, residuals)
        return residuals

    def set_dist_matrix (self, numseq, file):
        '''
        Sets the distance matrix for the "numseq" sequence.
        '''
        aln = AlignIO.read(open(file), 'stockholm')
        calculator = DistanceCalculator('blosum62')

        dist_matrix = calculator.get_distance(aln)
        matrix_list = self.__matrix_to_list(dist_matrix)

        if numseq == 1:
            self.matrix1 = matrix_list
            if self.matrix1.count(0.0) == len(self.matrix1):
                self.matrix1 = "zero"
        elif numseq == 2:
            self.matrix2 = matrix_list
            if self.matrix2.count(0.0) == len(self.matrix2):
                self.matrix2 = "zero"
        else:
            raise Exception ("numseq must be either 1 or 2!")


    def partial_corr(self):
        '''
        Computes partial correlation between "matrix1" and "matrix2"
        partialing out "sp_matrix"
        '''
        if self.sp_matrix is not None:
            res1 = self.__linear_reg_residues(self.matrix1, self.sp_matrix)
            res2 = self.__linear_reg_residues(self.matrix2, self.sp_matrix)
            self.partial_r = np.corrcoef(res1,res2)[1,0]
            return self.partial_r
        else:
            return "NA"


    def get_corr (self):
        '''
        Computes and returns all the possible correlation coefficients for this Interaction.
        '''
        if not self.matrix1 or not self.matrix2:
            raise Exception ("Distance matrix not calculated")
        elif self.matrix1 == "zero":
            raise Exception ("seq1.id matches 100%  with his homologs, distance matrix can't be calculated!")
        elif self.matrix2 =="zero":
            raise Exception("seq2.id matches 100%  with his homologs, distance matrix can't be calculated!")
        elif self.correlation is None:
            self.correlation = np.corrcoef(self.matrix1,self.matrix2)[1,0]
            self.spearman    = st.spearmanr(self.matrix1,self.matrix2)[0]
            return self.correlation, self.spearman
        else:
            return self.correlation, self.spearman

class ProgramNotFound(Exception):
    '''
    Exception thrown when a program needed by CheeseCake is not installed.
    '''
    def __init__(self, program):
        self.program = program

    def __str__(self):
        return "The program %s was not found in your system" %self.program
