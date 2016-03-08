from Bio import SeqIO
from Bio.Seq import Seq

class SeqRecordOrg(SeqIO.SeqRecord):
    '''
    This is a SeqRecord child class with species information
    '''
    def __init__(self, identifier, seq, name, species, description, dbxrefs, features, annotations, letter_annotations):
        super(SeqRecordOrg, self).__init__(
            id = identifier,
            seq = seq,
            name = name,
            description = description,
            dbxrefs = dbxrefs,
            features = features,
            annotations = annotations,
            letter_annotations = letter_annotations
        )
        self.species = species
