'''
Created on 04/02/2015

@author: flavio
'''
from Bio import Phylo
import networkx, pylab

apaf = Phylo.read('apaf.xml', 'phyloxml')

net = Phylo.to_networkx(apaf)
networkx.draw(net)
pylab.show()