'''
Created on 12/02/2015

@author: Flavio Lichtenstein
@local: UNIFESP DIS

http://pycogent.org/cookbook/using_likelihood_to_perform_evolutionary_analyses.html#reconstructing-ancestral-sequences

'''
from cogent import LoadTree, LoadSeqs
from cogent.evolve.models import F81
# import cogent.core.alignment.Alignment

tree = LoadTree('data/primate_brca1.tree')
aln = LoadSeqs('data/primate_brca1.fasta')
sm = F81()
lf = sm.makeLikelihoodFunction(tree, digits=3, space=2)
lf.setAlignment(aln)
lf.optimise(show_progress=False, local=True)

ancestors = lf.likelyAncestralSeqs()
print '\n--- nodes --------------------'
print ancestors

# for i in range(len(ancestors.AlignedSeqs)):
#    print "name: %s, seq: %s" % (ancestors.AlignedSeqs[i], ancestors.Names[i])

print ancestors.AlignedSeqs
print ancestors.Names


for key in ancestors.AlignedSeqs:
    print "name: %s, seq: %s" % (key, ancestors.AlignedSeqs[key]) 
n = 10
print '\n---- first %i lines of ancestral probabilities -------------'%(n)    

ancestral_probs = lf.reconstructAncestralSeqs()
print ancestral_probs['root'][:10]
print 'this is the root:', ancestors.AlignedSeqs['root']
