'''
Created on 12/02/2015

@author: http://pycogent.org/cookbook/using_likelihood_to_perform_evolutionary_analyses.html

@local: UNIFESP DIS

Evolutionary analysis using likelihood
Specifying substitution models
'''
# from gi.overrides.keysyms import lf

from cogent.evolve import models
print '\n-- model nucleotides '
print models.nucleotide_models
#['JC69', 'K80', 'F81', 'HKY85', 'TN93', 'GTR']
print '\n-- model codon '
print models.codon_models
#['CNFGTR', 'CNFHKY', 'MG94HKY', 'MG94GTR', 'GY94', 'H04G', 'H04GK', 'H04GGK']
print '\n-- model protein '
print models.protein_models
# ['DSO78', 'AH96', 'AH96_mtmammals', 'JTT92', 'WG01']

from cogent.evolve.models import F81
sub_mod = F81()

from cogent.evolve.models import GTR
sub_mod = GTR(with_rate=True, distribution='gamma')
print '\n\n-- model GTR '
print sub_mod

print '\n\n--- Making a likelihood function ----- '
from cogent import LoadTree
sub_mod = F81()
tree = LoadTree(treestring='(a,b,(c,d))')
lf = sub_mod.makeLikelihoodFunction(tree)

print '\n\n---  Providing an alignment to a likelihood function ----'
from cogent import LoadSeqs

sub_mod = F81()
tree = LoadTree(treestring='(a,b,(c,d))')
lf = sub_mod.makeLikelihoodFunction(tree)
aln = LoadSeqs(data=[('a', 'ACGT'), ('b', 'AC-T'), ('c', 'ACGT'),  ('d', 'AC-T')])
lf.setAlignment(aln)
print lf


print '\n\n--- Scoping parameters on trees ---- '
from cogent.evolve.models import CNFGTR
tree = LoadTree('data/primate_brca1.tree')
print tree.asciiArt()
sm = CNFGTR()
lf = sm.makeLikelihoodFunction(tree, digits=2)
lf.setParamRule('omega', tip_names=['Human', 'Orangutan'], outgroup_name='Galago', is_clade=True, init=0.5)

print lf


''' and much more '''
