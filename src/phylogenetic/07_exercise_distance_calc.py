'''
Created on 04/02/2015

@author: flavio
'''
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
# from Bio.Phylo.TreeConstruction import *


'''
In the csh shell: type 
setenv PATH "$PATH:/usr/local/bin/python" and press Enter.

In the bash shell (Linux): type 
export PATH="$PATH:/usr/local/bin/python" and press Enter.

In the sh or ksh shell: type 
PATH="$PATH:/usr/local/bin/python" and press Enter

import sys

for p in sys.path:
    print p 


'''


# aln = AlignIO.read('/usr/local/lib/python2.7/dist-packages/biopython-1.65-py2.7-linux-x86_64.egg/Bio/Tests/TreeConstruction/msa.phy', 'phylip')
aln = AlignIO.read('msa.phy', 'phylip')
print aln

calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)
print dm

constructor = DistanceTreeConstructor(calculator, 'nj')

tree = constructor.build_tree(aln)

print tree
''' 
supported_formats = {
    'newick': NewickIO,
    'nexus': NexusIO,
    'phyloxml': PhyloXMLIO,
    'nexml': NeXMLIO,
}'''
Phylo.write([tree], 'msa_nj.nhx', 'newick')

# print tree

tree.ladderize()   # Flip branches so deeper clades are displayed at top
draw = False
    

print '\n--- Leafs -------------'
leafs = tree.get_terminals()
for leaf in leafs:
    print leaf


print '\n---Common Ancestor -------------'
lca = tree.common_ancestor(leafs)
print lca
"""Most recent common ancestor (clade) of all the given targets.

Edge cases:
- If no target is given, returns self.root
- If 1 target is given, returns the target
- If any target is not found in this tree, raises a ValueError
"""


Phylo.draw(tree)
