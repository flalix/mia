'''
Created on 12/02/2015

@author: Flavio Lichtenstein
@local: UNIFESP DIS
'''

from cogent import LoadSeqs, DNA
from cogent.phylo import distance, nj
from cogent.draw.dendrogram import UnrootedDendrogram
from cogent.evolve.models import HKY85
import os

organism = 'Drosophila'
speciesList = ['americana', 'americana_americana', 'americana_texana', 'ananassae',  'arizonae', 'melanogaster', 'pseudoobscura', 'paulistorum', 'willistoni']
mnemList = ['ame', 'ameame', 'ametex', 'anana',  'arizo', 'mela', 'pseud', 'pauli'  ,'willi']


home = os.path.expanduser("~")
home = home.replace("\\","/")
root_files = home + '/data/%s/' %(organism)       
rootFasta = home + '/data/%s/fasta/'%(organism)
rootImage = home + '/data/%s/image/'%(organism)
rootTable = home + '/data/%s/files/'%(organism)
rootTree = home + '/data/%s/trees/'%(organism)

aln = None
i = -1

for spec in speciesList:
    i += 1
    mnem = mnemList[i]
    filename = '%s_maxmer_%s_Gene_Adh_100L_cutoff10_consensus'%(organism, spec)
    
    print filename
    
    if not aln:
        aln = LoadSeqs(rootFasta+"/"+filename+".fasta", moltype=DNA)
        key = aln.AlignedSeqs.keys()[0]
        try:
            keySimp = key[: key.find('_')]
        except:
            keySimp = key
            
        keySimp = mnem + '_' + keySimp
            
        aln = LoadSeqs(data= [(keySimp, aln.AlignedSeqs[key])], moltype=DNA)
        print len(aln.AlignedSeqs)
        
    else:
        aln2 = LoadSeqs(rootFasta+filename+".fasta", moltype=DNA)
        key = aln2.AlignedSeqs.keys()[0]
        try:
            keySimp = key[: key.find('_')]
        except:
            keySimp = key
            
        keySimp = mnem + '_' + keySimp

        aln2 = LoadSeqs(data= [(keySimp, aln2.AlignedSeqs[key])], moltype=DNA)
        try:
            aln = aln.addSeqs(aln2)
            print len(aln.AlignedSeqs)
        except:
            print 'error while merging'
            exit()
        
    
d = distance.EstimateDistances(aln, submodel= HKY85())
d.run()

mytree = nj.nj(d.getPairwiseDistances())
print mytree

filename = '%s_maxmer_Gene_Adh_100L_cutoff10_consensus'%(organism)
mytree.writeToFile(rootTree + filename + '.tree')

dendrogram = UnrootedDendrogram(mytree)
dendrogram.showFigure()


