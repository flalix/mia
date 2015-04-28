'''
Created on 19/02/2015

@author: Flavio
'''
import numpy as np


class MyEntropy(object):
    def __init__(self):
        pass
        
        
    def calcHSannon_bias_correction(self, seqs, numOfLetters):
               
        self.dicPiList = []
        self.HShannonList = []
        self.HShannonCorrList = []
        self.SeHShannonList = []
        self.SeHShannonCorrList = []

        lenSeq = len(seqs[0])
        numOfSeqs = len(seqs)
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
            hShanCorr = hShan + ((B-1)/ (2*N))
            
            VarCorr = 0
            ''' Roulston Eq. 40 '''
            for key in dicPi.keys():
                pi = dicPi[key]
                if (pi != 0.) and (pi != 1.):
                    VarCorr += (np.log(pi)+hShan)**2 * pi * (1-pi)
                
            SeHShannon     = hShan * np.sqrt(varPi/N)
            SeHShannonCorr = np.sqrt(VarCorr/N)
            
            self.dicPiList.append(dicPi)
            self.HShannonList.append(hShan)
            self.HShannonCorrList.append(hShanCorr)
            self.SeHShannonList.append(SeHShannon)
            self.SeHShannonCorrList.append(SeHShannonCorr)
                                    
    def crossCorrelation_ij(self, seqs, i,j, numOfLetters):
        dicPij = {}
        seqsLen = len(seqs)
        N = float(seqsLen)
        
        ''' scan all AGTC or AA 
            Var = 1/n * sum( pij *(1-pij) )

            SUM OF 3 FRAMES '''
        
        for row in range(seqsLen):
            n1 = seqs[row][i:i+numOfLetters] 
            n2 = seqs[row][j:j+numOfLetters] 

            try:                      
                dicPij[n1+n2] += 1
            except:
                dicPij[n1+n2] = 1
            
        for key in dicPij.keys():
            dicPij[key] = dicPij[key] / N  

        '''
            have been calculated at self.calcHSannon(seqs, numOfLetters)
            self.dicPiList.append(dicPi)
        '''
        dicPi = self.dicPiList[i]
        Bi = len(dicPi.keys())
             
        dicQj = self.dicPiList[j]
        Bj = len(dicQj.keys())
        
        Bij = len(dicPij.keys())
        
        correction = (Bi+Bj-Bij-1) / (2.*N)
        # print 'Bij %i,  Bi %i  Bj %i, and N = %i, correction=%f'%(Bij, Bi, Bj, N, correction)
        
        MItot = 0.
        varPij = 0.
        
        totVar = 0.
        totVarCorr = 0.
        hShanIJ = 0

        for key in dicPij.keys():
            n1 = key[0:numOfLetters]
            n2 = key[numOfLetters:numOfLetters+numOfLetters]
    
            Pij = dicPij[key]

            if Pij > 0 and Pij < 1:
                Pi = dicPi[n1]
                Qj = dicQj[n2]
                mi = Pij * np.log(Pij / (Pi * Qj) )
                    
                hShanIJ -= Pij * np.log(Pij)
                varPij = (1. + np.log(Pij))**2  * Pij * (1.-Pij)                
                totVar += varPij
                
                ''' Roulston equation 42 se = sqrt(var/N)'''
                varPijCorr =  (np.log(Pi) + np.log(Qj) - np.log(Pij) + mi)**2 * (Pij*(1.-Pij))
                totVarCorr += varPijCorr
                MItot += mi
                    

        if MItot < 0:
            print("*** Impossible *** MIk = %f, is less then 0, for i=%i  j=%i"%(MItot, i, j) )
                # raw_input('Enter to continue.')
        
        totVar = totVar*(hShanIJ**2) + N*(self.SeHShannonList[i]**2 + self.SeHShannonList[j]**2)
                            
        return MItot, MItot + correction, MItot*np.sqrt(totVar/N), np.sqrt(totVarCorr/N)
   
   
if __name__ == '__main__':
    
    ent = MyEntropy()
    
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
    seqs = seqs * 5
    L = len(seqs[0])
    
    '''
            self.dicPiList.append(dicPi)
            self.HShannonList.append(hShan)
            self.HShannonCorrList.append(hShanCorr)
            self.SeHShannonList.append(SeHShannon)
            self.SeHShannonCorrList.append(SeHShannonCorr)
    '''
                
    ent.calcHSannon_bias_correction(seqs, numOfLetters)


    print '--- dicPiList -----'
    for elem in ent.dicPiList:
        print elem
    
    print '\n--- HShannon -----'
    print np.round(ent.HShannonList,6)
    
    print '--- SE(HShannon) -----'
    print np.round(ent.SeHShannonList,6)
    
    print '\n--- HShannon Bias Correction -----'
    print np.round(ent.HShannonCorrList,6)
    
    print '--- SE(HShannon) Correction -----'
    print np.round(ent.SeHShannonCorrList,6)

    print '\n\n-------------'
    for i in range(0, L-1):
        for j in range(i+1, L):
            MI, MICorr, SeMi, SeMiCorr = \
            ent.crossCorrelation_ij(seqs, i,j, numOfLetters)
            
            print "MI=%f, SeMi=%f, MICorr=%f, SeMiCorr=%f" % (MI, SeMi, MICorr, SeMiCorr)

            