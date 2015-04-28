'''
Created on 19/02/2015

@author: Flavio
'''
import numpy as np

def calcHSannon(seqs, numOfLetters):
    dicPiList = []
    entropyShannon = []

    lenSeq = len(seqs[0])
    numOfSeqs = len(seqs)
    maxL = lenSeq-numOfLetters+1
    
    ''' for each nucleotide position (residue) '''
    for i in range(maxL):
        hShan = 0
        dicPi = {}
        
        ''' for each row '''
        for row in range(numOfSeqs):
            try:
                dicPi[seqs[row][i:i+numOfLetters]] += 1
            except:
                dicPi[seqs[row][i:i+numOfLetters]] = 1
            
        for key in dicPi.keys():
            dicPi[key] /= float(numOfSeqs) 
            ''' varPi += pi * (1-pi) * (np.log(pi) + 1.) ** 2 '''
            
            pi = dicPi[key]
            if (pi != 0.) and (pi != 1.):
                hShan -= pi * np.log(pi)
                
                
        dicPiList.append(dicPi)
        entropyShannon.append(hShan)
        
    return dicPiList, entropyShannon
    

seqs = []
numOfLetters = 1

seqs.append('AGTCA') 
seqs.append('ACTCC') 
seqs.append('AGCGG') 
seqs.append('AGTGT') 


dicPiList, entropyShannon =  calcHSannon(seqs, numOfLetters)

print '--- dicPiList -----'
print dicPiList

print '\n--- HShannon -----'
print entropyShannon


