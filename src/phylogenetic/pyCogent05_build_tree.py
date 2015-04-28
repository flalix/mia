'''
Created on 12/02/2015

@author: Flavio Lichtenstein
@local: UNIFESP DIS
'''

from cogent.draw.dendrogram import UnrootedDendrogram
from cogent import LoadSeqs, DNA
from cogent.app.fasttree import build_tree_from_alignment

''' load a tree
from cogent import LoadTree
tree = LoadTree('data/test.tree')
'''

aln = LoadSeqs('data/primate_brca1.fasta')

tree = build_tree_from_alignment(aln,DNA)

dendrogram = UnrootedDendrogram(tree)
dendrogram.showFigure()