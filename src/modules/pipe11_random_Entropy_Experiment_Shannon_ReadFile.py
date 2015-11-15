'''
Revised on 20/02/2014

@author: Flavio Lichtenstein
@local: UNIFESP - Bioinformatica
'''
import classes.BioPythonClass as bioClass

basic = bioClass.Basic()
showmessage = False

rootTable = 'C:/drosophila/files/'
randomShannonTableName = 'shannon_random_DNA_dic.txt'

entFile = bioClass.Entropy_File()

if not entFile.setNames(rootTable, randomShannonTableName):
    print 'Could not open ' + rootTable + randomShannonTableName
    
ok, msg = entFile.read_shannon_random_file(showKeys = False)

if not ok:
    print msg + ', perhaps Shannon dictionary bad formated.'
    exit()
   
    
print '--------------------------------------------------------'
print 'Random h-Shannon table opened.'
print '--------------------------------------------------------'
print 'L \t entropy \t sdv \t mer%'
print '--------------------------------------------------------'

for i in range(4, 500):
    ''' for DNA is i-3, 1st, 2nd .... '''
    posi = i-3
    randShannonEntropy, randShannonSdv = entFile.pos_shannon_random_file(posi, showmessage=False)
    '''
    print randShannonEntropy, randShannonSdv 
    print type(randShannonEntropy), type(randShannonSdv) 
    '''
    print i,'\t', round(randShannonEntropy,4), '\t', round(randShannonSdv,4),'\t', round( (randShannonEntropy/2)*100., 2)
print '--------------------------------------------------------\n\n'


