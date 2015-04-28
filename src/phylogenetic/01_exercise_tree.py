'''
Created on 04/02/2015

@author: flavio
'''

from Bio import Phylo


trees = Phylo.parse('phyloxml_examples.xml', 'phyloxml')

for tree in trees:
    print tree.name