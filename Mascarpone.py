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

class Distance_Calculator_Improved(DistanceCalculator):
    '''
    This distance calculator allows more letters (Y, N etc. for DNA, RNA etc.)
    '''
    def __init__(self,model="identity"):
        super(Distance_Calculator_Improved, self).__init__(model)

    def _pairwise(self, seq1, seq2):
        '''
        Calculate pairwise distance from two sequences.

        Returns a value between 0 (identical sequences) and 1 (completely
        different, or seq1 is an empty string.)
        '''
        score = 0
        max_score = 0
        if self.scoring_matrix:
            max_score1 = 0
            max_score2 = 0
            skip_letters = ['-', '*', '.', 'W', 'S', 'N', 'M', 'K', 'R', 'Y', 'B', 'V', 'D', 'H']
            for i in range(0, len(seq1)):
                l1 = seq1[i]
                l2 = seq2[i]
                if l1 in skip_letters or l2 in skip_letters:
                    continue
                if l1 not in self.scoring_matrix.names:
                    raise ValueError("Bad alphabet '%s' in sequence '%s' at position '%s'"
                                     % (l1, seq1.id, i))
                if l2 not in self.scoring_matrix.names:
                    raise ValueError("Bad alphabet '%s' in sequence '%s' at position '%s'"
                                     % (l2, seq2.id, i))
                max_score1 += self.scoring_matrix[l1, l1]
                max_score2 += self.scoring_matrix[l2, l2]
                score += self.scoring_matrix[l1, l2]
            # Take the higher score if the matrix is asymmetrical
            max_score = max(max_score1, max_score2)
        else:
            # Score by character identity, not skipping any special letters
            for i in range(0, len(seq1)):
                l1 = seq1[i]
                l2 = seq2[i]
                if l1 == l2:
                    score += 1
            max_score = len(seq1)

        if max_score == 0:
            return 1  # max possible scaled distance
        return 1 - (score * 1.0 / max_score)



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
        calculator = None
        if numseq == 3:
            calculator = Distance_Calculator_Improved('trans')
        else:
            calculator = DistanceCalculator('blosum62')

        dist_matrix = calculator.get_distance(aln)

        matrix_list = self.__matrix_to_list(dist_matrix)

        if numseq == 1:
            self.matrix1 = matrix_list
            print(len(self.matrix1))
            if self.matrix1.count(0.0) == len(self.matrix1):
                self.matrix1 = "zero"
        elif numseq == 2:
            self.matrix2 = matrix_list
            print(len(self.matrix2))
            if self.matrix2.count(0.0) == len(self.matrix2):
                self.matrix2 = "zero"
        elif numseq == 3:
            self.sp_matrix = matrix_list
            #print(dist_matrix)
            print (len(self.sp_matrix))
        else:
            raise Exception ("numseq must be either 1 or 2!")


    def partial_corr(self):
        '''
        Computes partial correlation between "matrix1" and "matrix2"
        partialing out "sp_matrix"
        '''
        if self.sp_matrix is not None:
            res1    = self.__linear_reg_residues(self.matrix1, self.sp_matrix)
            res2    = self.__linear_reg_residues(self.matrix2, self.sp_matrix)
            partial = np.corrcoef(res1,res2)[1,0]
            return partial
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
            self.partial_r   = self.partial_corr()
            return self.correlation, self.spearman, self.partial_r
        else:
            return self.correlation, self.spearman, self.partial_r

class ProgramNotFound(Exception):
    '''
    Exception thrown when a program needed by CheeseCake is not installed.
    '''
    def __init__(self, program):
        self.program = program

    def __str__(self):
        return "The program %s was not found in your system" %self.program
