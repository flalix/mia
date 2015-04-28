'''
Created on 04/02/2015

@author: flavio
'''
from Bio import Phylo

tree1 = Phylo.read('phyloxml01.xml', 'phyloxml')
print tree1
print '-------------'
tree2 = Phylo.read('phyloxml02.xml', 'phyloxml')
print tree2
print '-------------\n\n----Merging ---'

Phylo.write([tree1, tree2], 'phylo_merge0102.xml', 'phyloxml')

trees = Phylo.parse('phylo_merge0102.xml', 'phyloxml')

for t in trees:
    print t
    print '-----\n'


print '\nconverting from xml to newick'
print '----------------------------'

Phylo.convert('phylo_merge0102.xml', 'phyloxml', 'phylo_merge0102.nhx', 'newick')


trees = Phylo.parse('phylo_merge0102.nhx', 'newick')
print '\n---- Phylo newick tree'
for t in trees:
    print t
    print '-----\n'


print '\nconverting from xml to nexus'
print '----------------------------'

Phylo.convert('phylo_merge0102.nhx', 'newick', 'phylo_merge0102.nex', 'nexus')


trees = Phylo.parse('phylo_merge0102.nex', 'nexus')
print '\n---- Phylo nexus tree'
for t in trees:
    print t
    print '-----\n'




