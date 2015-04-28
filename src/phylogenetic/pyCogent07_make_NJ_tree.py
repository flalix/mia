'''
Created on 12/02/2015

@author: Flavio Lichtenstein
@local: UNIFESP DIS
'''

from cogent import LoadSeqs
from cogent.phylo import distance, nj
from cogent.draw.dendrogram import UnrootedDendrogram
from cogent.evolve.models import HKY85

al = LoadSeqs("data/long_testseqs.fasta")

d = distance.EstimateDistances(al, submodel= HKY85())
d.run()

mytree = nj.nj(d.getPairwiseDistances())
print mytree

mytree.writeToFile('data/test_nj.tree')


dendrogram = UnrootedDendrogram(mytree)
dendrogram.showFigure()

