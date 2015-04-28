'''
Created on 27/11/2014

@author: Flavio Lichtenstein
'''
import numpy as np
import classes.Timing as timePack
from tests import test_Class as TC

numOfLetters = 2
isProtein = False

tc = TC.test_Class(isProtein, numOfLetters) 

def calcVerticalPi(seqs):
    Pis = []
    entropyShannon = []

    lenSeq = len(seqs[0])
    numOfSeqs = len(seqs)
    max_len = lenSeq-numOfLetters+1
    
    print '--------'
    print tc.nuc_aa_List
    print '--------'
    
    for i in range(max_len):
        listTot = []
        hShan = 0
        
        dicPi = {}
        for n1 in tc.nuc_aa_List:
            dicPi[n1] = 0
        
        for row in range(numOfSeqs):
            dicPi[seqs[row][i:i+numOfLetters]] += 1
            
        for n1 in tc.nuc_aa_List:
            listTot.append(dicPi[n1])
            dicPi[n1] /= float(numOfSeqs) 

            pi = dicPi[n1]
            if (pi != 0.) and (pi != 1.):
                hShan -= pi * np.log(pi)
                
                
        Pis.append([listTot, dicPi])
        entropyShannon.append(hShan)
        
    return Pis, entropyShannon

        
time = timePack.Timing()    
seqs = []
for _ in range(10000):
    seqs.append('AAAGTC')        
for _ in range(2000):
    seqs.append('AAACTC')        
for _ in range(3000):
    seqs.append('AAAGAC')        

_ = time.start()    
Pis, hShan = calcVerticalPi(seqs)
_ = time.finish()
print time.milli(), 'ms\n\n'

'''
for elem in Pis:
    print '\t', elem[0]
    print '\t', elem[1]
print '------------------\n'
'''
   
for i in range(len(Pis)):
    elem = Pis[i]
    print 'from %i to %i'%(i,i+numOfLetters-1)
    print '\t', elem[0]
    print '\t', elem[1]
print '------------------\n'
print hShan


