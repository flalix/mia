'''
Created on Mar 16, 2015

@author: flavio
'''
import classes.Timing as timePack

def cut_vertical_gaps(cut, seqs_align):
    cutoff = cut/100.
    
    seqs = []
    for seq in seqs_align:
        seqs.append(seq)
        
    seqLen = len(seqs[0])
    numSeqs = len(seqs)
    
    minGap = int(cutoff * numSeqs)
    gaps = []
    
    for col in range(seqLen):
        count = 0
        for row in range(numSeqs):
            if seqs[row][col] <> '-':
                count += 1
                
        if count < minGap:
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

        
    for row in range(numSeqs):
        seq = seqs[row]
        seq2 = ''
        for limInf, limSup in replaces:
            seq2 += seq[limInf: limSup]
            
        seqs_align[row] = seq2
    
    return seqs_align

    
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

        