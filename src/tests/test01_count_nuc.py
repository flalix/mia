'''
Created on 24/11/2014

@author: Flavio Lichtenstein
'''
import classes.Timing as timePack
from tests import test_Class as TC

tc = TC.test_Class()

time = timePack.Timing()    

striNuc = 'AGTCAAGG'*1000
seqs = []
for i in range(100):
    seqs.append(striNuc)

_ = time.start()
ni, pi = tc.calc_DNA_pi(seqs)
_ = time.finish()
print time.milli(), 'ms'

print 'ni=%s  pi=%s' %(str(ni), str(pi))


#---- amino acids ----------

striAA = 'AGCLLLWLWL'*1000
seqs = []
for i in range(100):
    seqs.append(striAA)


_ = time.start()
ni, pi = tc.calc_Prot_pi(seqs)
_ = time.finish()
print time.milli(), 'ms'

print 'ni=%s  pi=%s' %(str(ni), str(pi))


