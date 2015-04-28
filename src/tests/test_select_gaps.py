'''
Created on 26/01/2015

@author: flavio
'''



seqLen = 403
gaps = [0,1,2,3,7, 8,9,13, 400,401]

seqLen = 401
gaps = [0,1,2,3,7, 8,9,13, 400,401]

seqLen = 430
gaps = [2,3,7, 8,9,13, 400,401]


this = -1
replaces = []

for k in gaps:
    this += 1
    
    if k == this:
        continue

    replaces.append([this,k])
    this = k

this = k + 1
if this < seqLen:
    replaces.append([this,seqLen])
    

print replaces
    