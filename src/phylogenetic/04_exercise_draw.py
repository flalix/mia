'''
Created on 04/02/2015

@author: flavio
'''
from Bio import Phylo
import pylab


tree = Phylo.read('phyloxml01.xml', 'phyloxml')
print tree

num = 3

if num == 1:
    tree.ladderize()   # Flip branches so deeper clades are displayed at top
    Phylo.draw(tree)
elif num == 2:
    Phylo.draw_ascii(tree)
else:    
    # sudo pip install networkx
    print 'graphviz ...'
    Phylo.draw_graphviz(tree)
    pylab.show()
