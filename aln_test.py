import argparse
import os
import sys
import socket
import time
import numpy as np
from Bio import SeqIO
from Bio import AlignIO
from Bio.Blast import NCBIXML
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
from Bio.Align.Applications import TCoffeeCommandline     as tcoffee

#-----------------------------------------------------------------------

aln = AlignIO.read(open('tmp/1.fa.aln'), 'clustal')
calculator = DistanceCalculator('blosum62')
dist_matrix = calculator.get_distance(aln)

i=0
j=0
da_list = list()
for row in dist_matrix:
    print ('New Row!')
    j=0
    for column in row:
        if i<j:          # with this, you take out the 0's so n = (N²-N)/2
            print (dist_matrix[i,j])
            da_list.append(dist_matrix[i,j])
            
        j+=1
    i+=1
print (da_list)
print(len(da_list))
result = np.corrcoef(da_list,da_list)[1,0]
print("Result: "+str(result))



