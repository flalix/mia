'''
Created on 24/11/2014

@author: Flavio Lichtenstein
'''
import classes.Timing as timePack
import classes.BioPythonClass as bioClass
from tests import test_Class as TC
import classes.Drosophila as dro
import numpy as np

tc = TC.test_Class()
time = timePack.Timing()    
dr = dro.Drosophila()

class Desk:
    def __init__(self):
        self.basic = bioClass.Basic()

        self.withCorrection = True
        self.showmsg_obs = None
        
        self.numOfLetters = 1
        self.organism = 'Drosophila'
        self.seqType = 'Gene'
        self.gene = 'Adh'
        self.frame = None
        self.cutoffLength = 400
        self.cutoffNumSeq = 10
        self.isProtein = False

        self.species = 'melanogaster'
        
        self.filename = ('%s_%s_%s_%s_%s_%iL_cutoff%i_consensus.fasta') %\
             (self.organism, 'mincut', self.species, self.seqType, self.gene, \
              self.cutoffLength, self.cutoffNumSeq)
        
        self.name = ('%s_%s') %  (self.organism, self.species)
        self.names = [[self.filename, self.name]]

        self.de_novo = True
        self.showmessage = False
        
        self.all_species = 'mincut'
        
        self.rootFasta = 'C:/data/fasta2014/'
        self.rootImage = 'C:/data/image/'
        self.rootTable = 'C:/data/files/'
        
desk = Desk()
desk.nuc_aa_List = tc.nuc_aa_List

lista2 = [desk.species]
vmi = TC.Vertical_MI(desk, maxFig=1, showmessage=False)
            
if not vmi.ok:
    print 'Problems with Vertical_MI.'
    exit()
 
if not vmi.read_file_new(desk.names[0]):
    print "Could not find %s"%(desk.names[0])
    exit()

try:
    numOfSeqs = len(vmi.ent.mySequence.sequence)
except:
    print "error finding sequences in %s"%(desk.filename)
    exit()

iLoop = 0

stri = '%i/%i) %s - %s and aligned'%(iLoop, len(lista2), desk.species, 'mincut') 
   
    
if desk.organism == 'Drosophila':
    sp = dr.mnemonic(desk.species.replace('Drosophila ',''))
else:
    sp = desk.species
                
'''
if len(sp) < 12:
    sp = sp + ' '*(12-len(desk.species))
speciesList.append([desk.species,sp,numOfSeqs])
'''
   
seqs = []
for rec in vmi.ent.mySequence.sequence:
    seqs.append(str(rec)) 

_ = time.start()    
    
L = len(seqs[0])
nSeqs = len(seqs)
                        
'''  all params in mySequence '''
vmi.ent.mySequence.set_seq_all_param(desk, seqs)
stri +=  ' >> #seqs: %i, len=%i' % (numOfSeqs, len(seqs[0]) )
print stri

sufix = vmi.ent.mySequence.sufix_vmi
ret, msg = vmi.prepareVerticalMutualInfo(desk, iLoop, sufix=sufix, showmessage=desk.showmessage)
if not ret:
    print msg
    exit()

if not desk.withCorrection:           
    arrayMI = np.array(vmi.listMI)
    arraySE = np.array(vmi.listVar)
else:
    arrayMI = np.array(vmi.listMIinf)
    arraySE = np.array(vmi.listVar)
    
arraySE = np.sqrt( arraySE / nSeqs)

_ = time.finish()
print time.milli(), 'ms'




