'''
Created on Mar 16, 2015

@author: flavio
'''

'''
Created on Mar 16, 2015

@author: flavio
'''
import classes.Timing as timePack
import numpy as np

'''
a = ['a','b','a','c','a','a']
a = np.array(a)

a = np.array(a)
print list(a).count('a')
exit()
'''

def cut_vertical_gaps(cutoff, seqs_align):
    cutoff = cutoff/100.
    
    seqs = []
    for seq in seqs_align:
        seqs.append([x for x in seq])

    seqs = np.array(seqs)

    seqLen = len(seqs[0])
    numSeqs = len(seqs)
        
    maxGap = numSeqs - int(cutoff * numSeqs)
    gaps = []
    
    for col in range(seqLen):
        count = list(seqs[:,col]).count('-')

        if count > maxGap:
            gaps.append(col)

    if not gaps:
        return seqs_align
    
    this = -1
    replaces = []
    
    for k in gaps:
        this += 1
        
        if k == this:
            continue
    
        replaces.append([this,k])
        this = k
    
    this = gaps[len(gaps)-1] + 1
    if this < seqLen:
        replaces.append([this,seqLen])
        
    
    print 'Selected bp regions', replaces

    return [ (''.join( [seqs_align[row][limInf: limSup] for limInf, limSup in replaces] )) for row in range(numSeqs) ]
    '''
    for row in range(numSeqs):
        seq = seqs[row]
        seq2 = ''
        for limInf, limSup in replaces:
            seq2 += seq[limInf: limSup]
            
        seqs_align[row] = seq2
   
    return seqs_align
    '''
    
        
    
k = 10
seqs = []
seqs.append('1AGTCABCAGTCABCAGTCABC'*k)
seqs.append('2AGTCABCAGTCABCAGTCABC'*k)
seqs.append('3AGTCABCAGTCABCAGTCABC'*k)
seqs.append('4AGT--BCAGTCABCA--CABC'*k)
seqs.append('5AGT--BCAGTCABCA--CABC'*k)
seqs.append('6AGT--BCAGTCABCA--CABC'*k)
seqs.append('7AGT------TCAB----CABC'*k)
seqs.append('8AGTCABCAGTCABCAGTCABC'*k)
seqs.append('9AGT------TCAB----CABC'*k)
seqs.append('0AGTCABCAGTCABCAGTCABC'*k)

seqs *= 10000
print('len: %i', len(seqs))
print('L: %i', len(seqs[0]))

cutoff = 75
print("cutoff:", cutoff)

time = timePack.Timing()
time.start()
seqs = cut_vertical_gaps(cutoff, seqs)
time.finish()
print time.milli(), 'ms'

print('L: %i', len(seqs[0]))
