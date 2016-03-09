from Bio import SeqIO
from Bio.Seq import Seq

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
        self.homologs     = dict()
        self.__homolog_sp = None

    def get_homolog_species(self):
        '''
        This returns a set with the species of all the homologs of this progein
        '''
        if self.__homolog_sp is None:
            self.__homolog_sp = set()
            for homo_name, homo_seq in self.homologs.items():
                self.__homolog_sp.add(homo_seq.species)
            return self.__homolog_sp
        else:
            return self.__homolog_sp
