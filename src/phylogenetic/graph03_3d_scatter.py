'''
Created on 26/02/2015

@author: flavio
'''
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

'''
import math
a= int(math.ceil(4.01))
print a, type(a)
exit()

import math
lista = [5, 10, 15, 35, 40, 67, 83.6, 95, 100]
delta = 20.

for i in range(len(lista)):
    w = int(math.ceil(lista[i]/delta)) - 1
    
    print i, lista[i], w
    
exit()
'''

def randrange(n, vmin, vmax):
    return (vmax-vmin)*np.random.rand(n) + vmin

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
n = 100
for c, m, zl, zh, cutoff in [('r', 'o', 0, 20, 15), ('b', '^', 0, 40, 30)]:

    print 'random'
    xs = list(randrange(n, 20, 50)) + list(randrange(n, 120, 150))
    ys = list(randrange(n, 40, 100)) + list(randrange(n, 200, 220))
    zs = list(randrange(2*n, zl, zh))
    
    print len(zs)
    print 'calc low items'
    
    ''''    
    low = [i for i in range(len(zs)) if zs[i] < cutoff]
    
    print 'rebuild with list'
    x1, y1, z1 = [],[],[]
    for i in range(len(xs)):
        if i not in low:
            x1.append(xs[i])
            y1.append(ys[i])
            z1.append(zs[i])
    '''
        
    low = {i:None for i in range(len(zs)) if zs[i] < cutoff}
    
    print 'rebuild with list dic'

    x1, y1, z1 = [],[],[]
    for i in range(len(xs)):
        if i not in low:
            x1.append(xs[i])
            y1.append(ys[i])
            z1.append(zs[i])
                
    ax.scatter(x1, y1, z1, c=c, marker=m)
     

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()

'''
    print 'cleaning low items'
    for i in range(len(low)):
        xs.pop(low[i]-i)
        ys.pop(low[i]-i)
        zs.pop(low[i]-i)
    '''            