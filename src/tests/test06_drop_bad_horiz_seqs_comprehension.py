'''
Created on Mar 16, 2015

@author: flavio
'''
import classes.Timing as timePack

def cut_horizotal_gaps( minimum, seqs):
    cutoff = 1. - minimum/100.
    if cutoff == 0:
        return seqs

    L = len(seqs[0])
    numSeqs = len(seqs)
    
    maxGap = int(L * cutoff)

        
    badList = [ i for i in range(numSeqs) if  seqs[i].count('-') > maxGap]
               
    '''
    badList = []
    for i in range(numSeqs):
        count = seqs[i].count('-')
                
        if count > cutGap:
            badList.append(i)
    '''    

    print(maxGap, badList)
    offset = 0
    for k in badList:
        seqs.pop(k-offset)
        offset += 1
    
    return seqs


        
seqs = []
seqs.append('1AGTCABCAGTCABCAGTCABC')
seqs.append('2AGTCABCAGTCABCAGTCABC')
seqs.append('3AGTCABCAGTCABCAGTCABC')
seqs.append('4AGT--BCAGTCABCA--CABC')
seqs.append('5AGT--BCAGTCABCA--CABC')
seqs.append('6AGT--BCAGTCABCA--CABC')
seqs.append('7AGT------TCAB----CABC')
seqs.append('8AGTCABCAGTCABCAGTCABC')
seqs.append('9AGT------TCAB----CABC')
seqs.append('0AGTCABCAGTCABCAGTCABC')

seqs *= 1 
L = len(seqs[0])
print('L: %i'%L)
print('len: %i'%len(seqs))

cutoff = 75
print("cutoff: %i = %i"%(cutoff, int(L * (1-cutoff /100.))) )


if len(seqs) < 30:
    for seq in seqs:
        print seq
    
time = timePack.Timing()
time.start()

seqs = cut_horizotal_gaps(cutoff, seqs)

time.finish()
print time.milli(), 'ms'



if len(seqs) < 30:
    for seq in seqs:
        print seq

