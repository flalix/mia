'''
Created on 16/09/2015

@author: Flavio
'''
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


def shuffling(x):
    if type(x)==str:
        x = [a for a in x]
        random.shuffle(x)   
        x = ''.join(str(e) for e in x)   
        return x
    
    if type(x)==int:
        if x < 0:
            isNegative = True
            x = -x
        else:
            isNegative = False
            
        x = [a for a in str(x)]
        random.shuffle(x)   
        x = ''.join(str(e) for e in x) 
        return -int(x) if isNegative else int(x)
    
    if type(x)==float:
        if x < 0:
            isNegative = True
            x = -x
        else:
            isNegative = False
            
        x = [a for a in str(x)]
        random.shuffle(x)   
        x = ''.join(str(e) for e in x)   
        return -float(x) if isNegative else float(x)  
    
    if type(x)==list:
        random.shuffle(x)   
        return x

    if type(x)==Seq:
        seq = str(x)
        seq = [x for x in seq]
        random.shuffle(seq)   
        seq = ''.join(str(e) for e in seq)  
        x = Seq(seq, IUPAC.unambiguous_dna)
        return x

    if type(x)==SeqRecord:
        idi = x.id
        desc = x.description
        name = x.name

        seq = str(x.seq)
        seq = [x for x in seq]
        random.shuffle(seq)   
        seq = ''.join(str(e) for e in seq)  
        x = Seq(seq, IUPAC.unambiguous_dna)
        return SeqRecord(x, id=idi, description=desc, name=name)
    
    print("Could not shuffle")
    return x
    


print shuffling("abcde")
print shuffling(list(range(16)))
print shuffling(Seq('AACCGGTTCCAA', IUPAC.unambiguous_dna))
print shuffling(SeqRecord(Seq('AACCGGTTCCAA', IUPAC.unambiguous_dna),
                                id="any", description='') )
print shuffling(2342)
print shuffling(-456)
print shuffling(18.379)
print shuffling(-0.056)

'''
x = [i for i in range(30)]
print x
random.shuffle(x)  
print x
print type(x), type(x)==list
L = 50
seq = [random.sample("ACGT", 1)[0] for _ in range(L)]


seq = ''.join(str(e) for e in seq)
print(seq)  
seq = [x for x in seq]
random.shuffle(seq)   
seq = ''.join(str(e) for e in seq)   
print(seq)


rec = SeqRecord(Seq(seq, IUPAC.unambiguous_dna),
                                id="any", description='')
print "printing rec"
print rec
print rec.seq

print "----------------"
seq = str(rec.seq)
seq = [x for x in seq]
random.shuffle(seq)   
seq = ''.join(str(e) for e in seq)  

seq = Seq(seq, IUPAC.unambiguous_dna)

rec = SeqRecord(seq, id="any", description='')

print(rec.seq)
print type(seq), type(seq)==Seq
print type(rec), type(rec)==SeqRecord

'''


