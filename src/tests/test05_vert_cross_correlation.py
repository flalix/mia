'''
Created on 27/11/2014

@author: Flavio Lichtenstein
'''
import numpy as np
import classes.Timing as timePack
from tests import test_Class as TC

numOfLetters = 1
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

def crossCorrelation_ij(seqs, i,j):
    dicPij = {}
    dicCij = {}
    
    seqLen = len(seqs)
    ''' scan all AGTC or AA 
        Var = 1/n * sum( pij *(1-pij) )

        SUM OF 3 FRAMES '''
    
    for row in range(seqLen):
        n1 = seqs[row][i:i+numOfLetters] 
        n2 = seqs[row][j:j+numOfLetters] 

        try:                      
            dicCij[n1+n2] += 1
        except:
            dicCij[n1+n2] = 1
        
    for key in dicCij.keys():
        dicPij[key] = dicCij[key] / float(seqLen)  


    '''
        have been calculated at calcVerticalPi(seqs, nuc_aa_List, numOfLetters)
        Pis.append([listTot, dicPi])
    '''
    print i,j
    print '\t', Pis[i]
    print '\t', Pis[j]
    ciAux = np.array(Pis[i][0])
    dicPi = Pis[i][1]
    
    Bi = len(ciAux[ np.where( ciAux > 0) ])
         
    cjAux = np.array(Pis[j][0])
    dicQj = Pis[j][1]

    Bj = len(cjAux[ np.where( cjAux > 0) ])
    
    cij = np.array(dicCij.values())
    Bij = len(cij[ np.where( cij > 0) ])
    # totCount = np.sum(cij)

    correction = (Bi+Bj-Bij-1) / (2.* float(seqLen))
    # print 'Bij %i,  Bi %i  Bj %i, and N = %i, correction=%f'%(Bij, Bi, Bj, N, correction)
    
    totVar = 0
    '''
    print 'nuc1,\tnuc2,\tPij,\t\tPi,\t\tQj,\t\tmi,\t\tMIk'
    '''
    MItot = 0.
    for key in dicPij.keys():
        n1 = key[0]
        n2 = key[1]

        Pij = dicPij[key]

        print '\t\t ', key, ' pi',dicPi[n1],' Qj',dicQj[n2], ' Pij', Pij,

        if Pij > 0 and Pij < 1:
            Pi = dicPi[n1]
            Qj = dicQj[n2]
            mi = Pij * np.log(Pij / (Pi * Qj) )
            print 'mi', mi
            
            ''' Roulston equation 42 se = sqrt(var/N)'''
            var =  (np.log(Pi) + np.log(Qj) - np.log(Pij) + mi)**2 * (Pij*(1.-Pij))
            totVar += var
            MItot += mi
        else:
            print ''
                

    if MItot < 0:
        print("*** Impossible *** MIk = %f, is less then 0, for i=%i  j=%i"%(MItot, i, j) )
            # raw_input('Enter to continue.')
    
    print '\t', MItot, MItot + correction, totVar
    return MItot, MItot + correction, totVar

time = timePack.Timing()    
seqs = []
for _ in range(10000):
    seqs.append('AAAGTC')        
for _ in range(2000):
    seqs.append('AAACTC')        
for _ in range(3000):
    seqs.append('AAAGAC')        

Pis, hShan = calcVerticalPi(seqs)
max_len = len(seqs[0]) - (numOfLetters-1)

_ = time.start()

for i in range(max_len-1):
    for j in range(i+1, max_len ):
        MIk, MIkCorr, varMI = crossCorrelation_ij(seqs, i,j)
        
_ = time.finish()
print time.milli(), 'ms\n\n'

print i,j,MIk
