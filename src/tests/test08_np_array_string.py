'''
Created on Mar 16, 2015

@author: flavio
'''
import numpy as np

a = [['a','a','b','c','d','a','a'], ['a','a','a','a','a','a','a'], ['a','a','b','c','d','a','a']]
print a
a = np.array(a)

print a
print a[:,2:5]

a = 'ABCDEFGHIJKLM'*10
a = [x for x in a]
print(a)

a = [a,a,a,a]
a = np.array(a)
print a[:,2:8]

exit()
