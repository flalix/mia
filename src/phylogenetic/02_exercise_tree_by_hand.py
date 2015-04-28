'''
Created on 04/02/2015

@author: flavio
'''
from Bio import Phylo

from cStringIO import StringIO
 
treedata = "(A, (B, C), (D, E))"
handle = StringIO(treedata)
tree = Phylo.read(handle, "newick")

print tree
