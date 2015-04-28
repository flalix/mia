'''
Created on 04/02/2015

@author: flavio
'''
from Bio import Phylo


apaf = Phylo.read('apaf.xml', 'phyloxml')

Phylo.draw_ascii(apaf)
