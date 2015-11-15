'''
Created on 02/09/2013

@author: Flavio Lichtenstein
@local: UNIFESP - Bioinformatica
'''
from Bio import Entrez
import classes.Drosophila as cDro

Entrez.email = "flalix@gmail.com"

db="nucleotide"
showmessage = False
root = 'C:/drosophila/fasta2013/'

fileNameSearch = 'C:/drosophila/files/drosophila_search.dic'
fileNameLog = 'C:/drosophila/files/drosophila_search__log.dic'

'''
listOrSpecies = ''
for specie in drosGene.listSpecies:
    listOrSpecies += ' or ' + specie + '[Organism]'
'''

def save_log(dic):
    f = open(fileNameLog, 'w')
    
    s = str(dic)
    f.write(s)
    f.flush()
    f.close()

def save_data(dicSpecGene, cutoff):
    dic = {}
    
    for gene in dicSpecGene.keys():
        if len( dicSpecGene[gene]) > cutoff:
            dic[gene] = dicSpecGene[gene]
            
    f = open(fileNameSearch, 'w')
    
    s = str(dic)
    f.write(s)
    f.flush()
    f.close()
            
''' already done '''
try:
    f = open(fileNameSearch, 'r')
    dicSpecGene = eval(f.read())
    f.close()
    found = True
    
except:
    found = False

cutoff = 9
            
if not found:
    ''' log save partial written data '''
    try:
        f = open(fileNameLog, 'r')
        dicSpecGene = eval(f.read())
        f.close()
    
    except:
        dicSpecGene = {}
    
    genes = ['run'] # drosGene.listGenes
    
    count = 0
    for gene in genes:
        count += 1

        if gene in dicSpecGene.keys():
            print '*', count, '/', len(genes), gene, ':', dicSpecGene[gene]
            continue
        
        print count, '/', len(genes), gene, ':',
        
        dicSpecGene[gene] = []
        
        j = 0
        itsOk = True
        repeat = 0
        getOut = True
        
        dro = cDro.Drosophila()
        
        while (True):
            for Specie in dro.speciesList:
                try:
                    j += 1
                    print j,
                    searchTerm = 'Drosophila %s[Organism] %s[Gene] not shotgun[Title] not complete sequence[Title] not complete genome[Title] not chromosome[Title] not genomic contig[Title] not supercontig [Title]'%(Specie, gene)
                    #if count == 753:
                    #    print searchTerm
                    
                    handle = Entrez.esearch(db="Nucleotide", term=searchTerm) # , rettype='gb', retmax=100000)
                    record = Entrez.read(handle)
                    
                    num = int(record["Count"])  # len(record["IdList"])
                    #print num
                    
                    if num > cutoff:
                        print ' [', Specie, num,']',
        
                    dicSpecGene[gene].append([Specie, num])
                    getOut = True
                    itsOk = True
                
                except:
                    j = 0
                    if repeat == 0:
                        print 'again',
                    repeat += 1 
                    
                    if repeat == 20:
                        itsOk = False
                        getOut = True
                        break
                    else:
                        getOut = False
                        print '(', repeat+1, ')',
                        break

            
            if getOut:
                break
            
        if not itsOk:
            print '\n\n------------ problems exceed 20 times search ------------'
            exit()
            
        save_log(dicSpecGene)
        print ''

    save_data(dicSpecGene, cutoff)

if not found:
    try:
        f = open(fileNameSearch, 'r')
        dicSpecGene = eval(f.read())
        f.close()
        
    except:
        print 'problems reading', fileNameSearch
        exit()

goodDic = {}

for gene in dicSpecGene.keys():
    has3orMore = 0
    filterElem = []
    for elem in  dicSpecGene[gene]:
        if elem[1] > cutoff:
            has3orMore += 1
            filterElem.append(elem)

    if has3orMore > 2:
        goodDic[gene] = filterElem
         
        
print 'Gene', '\t', 'Species', '\t', 'Num of Seqs'
for gene in goodDic.keys():
    print gene, '\t', 
    for elem in  goodDic[gene]:
        print elem[0], '\t', elem[1], '\t',
    print ''        

'''
    filenameGBK = root+'drosophila_CDS_ONLY_%s.gbk' % (gene)
    print gene, filenameGBK
'''
print '\n\n ------------------------------------ \n'
lista = list(goodDic.keys())
lista.sort()

for gene in lista:
    print "'"+gene + "',",
    
        