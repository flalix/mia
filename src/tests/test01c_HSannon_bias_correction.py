'''
Created on 19/02/2015

@author: Flavio
'''
import numpy as np

def calcHSannon_bias_correction(seqs, numOfLetters):
    dicPiList = []
    entropyShannon = []
    entropyShannonCorr = []
    seEntropyShannon = []
    seEntropyShannonCorr = []

    lenSeq = len(seqs[0])
    numOfSeqs = len(seqs)
    print 'Begining, numOfSeqs=%i, lenSeq=%i\n'%(numOfSeqs,lenSeq)
    maxL = lenSeq-numOfLetters+1
    
    for i in range(maxL):
        dicPi = {}
        varPi = 0
        hShan = 0
        
        for row in range(numOfSeqs):
            try:
                dicPi[seqs[row][i:i+numOfLetters]] += 1
            except:
                dicPi[seqs[row][i:i+numOfLetters]] = 1
            
        for key in dicPi.keys():
            dicPi[key] /= float(numOfSeqs) 
            
            pi = dicPi[key]
            if (pi != 0.) and (pi != 1.):
                hShan -= pi * np.log(pi)
                varPi += (1. + np.log(pi))**2  * pi * (1.-pi)
                
        B = len(dicPi.keys())
        N = float(numOfSeqs)
        
        ''' Roulston Eq. 13 '''
        Hcorr = hShan + ((B-1)/ (2*N))
        
        VARinf = 0
        if i == 4:
            pass
        ''' Roulston Eq. 40 '''
        for key in dicPi.keys():
            pi = dicPi[key]
            if (pi != 0.) and (pi != 1.):
                VARinf += (np.log(pi)+hShan)**2 * pi * (1.-pi)
            
        SE     = hShan * np.sqrt(varPi/N)
        SEcorr = np.sqrt(VARinf/N)
        
        dicPiList.append(dicPi)
        entropyShannon.append(hShan)
        entropyShannonCorr.append(Hcorr)
        seEntropyShannon.append(SE)
        seEntropyShannonCorr.append(SEcorr)
        
    return dicPiList, entropyShannon, entropyShannonCorr, seEntropyShannon,seEntropyShannonCorr 

seqs = []
numOfLetters = 1

seqs.append('AGTCA') 
seqs.append('ACTCC') 
seqs.append('AGCGG') 
seqs.append('AGTGT') 

'''
for i in range(90000):
    seqs.append('AGTCA') 
'''
seqs = seqs * 300

dicPiList, entropyShannon, entropyShannonCorr, \
seEntropyShannon,seEntropyShannonCorr = \
   calcHSannon_bias_correction(seqs, numOfLetters)

print '--- dicPiList -----'
for elem in dicPiList:
    print elem

print '\n--- HShannon -----'
print np.round(entropyShannon,6)

print '--- SE(HShannon) -----'
print np.round(seEntropyShannon,6)

print '\n--- HShannon Bias Correction -----'
print np.round(entropyShannonCorr,6)

print '--- SE(HShannon) Correction -----'
print np.round(seEntropyShannonCorr,6)

