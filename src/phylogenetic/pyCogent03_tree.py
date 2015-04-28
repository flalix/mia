'''
Created on 12/02/2015

@author: Flavio Lichtenstein
@local: UNIFESP DIS
'''
from cogent import LoadSeqs
from cogent.phylo import distance
from cogent.evolve.models import F81

aln = LoadSeqs('data/primate_brca1.fasta')
d = distance.EstimateDistances(aln, submodel=F81())
d.run()

dist_opt_args = dict(max_restarts=5, max_evaluations=10000)
d.run(dist_opt_args=dist_opt_args)
print d

from cogent.phylo import nj
njtree = nj.nj(d.getPairwiseDistances())
njtree = njtree.balanced()
print njtree.asciiArt()

