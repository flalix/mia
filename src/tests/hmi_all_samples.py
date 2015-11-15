'''
Created on Sep 14, 2015

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica
'''
import tests.Desktop_test as Desktop_test
import random
import classes.Drosophila as dro

organism = "Drosophila"      
gene = "Adh"      
title = ''
cutoffLength = 100
cutoffNumSeq = 7

filename_default = 'default.ini'
        
if organism == "Drosophila":
    dr = dro.Drosophila()
else:
    dr = None
    
desk = Desktop_test.Desktop(filename_default, organism, gene, title, cutoffLength, cutoffNumSeq)

opt = "mincut"
desk.minmax =  opt

desk.NtreshSamples = 10
desk.NminSamples = 6
desk.sampleElems = 4

desk.de_novo = False
desk.showGraph = False
desk.saveGraph = False
desk.saveData = False

#----- fix here ------------
desk.simulation = False
#---------------------------
shuffle = False
if shuffle: simulattion = False

print desk.filename_species_list

desk.build_listSpecies()
print desk.cutoffNumSeq, desk.totList

desk.dic_sel = {}
desk.dic_not_sel = {}

for key in desk.dicParams.keys():
    mat = desk.dicParams[key]
    # print key, mat

    Nmincut = mat[2]
    Nmaxmer = mat[3]
    if opt == "mincut":
        N = Nmincut
    else:
        N = Nmaxmer 
    
    if N >= desk.cutoffNumSeq:
        desk.dic_sel[key] = mat
    elif N >= desk.NminSamples:
        desk.dic_not_sel[key] = mat
        

selSpecies = {}

Lmincut = -1
Lmaxmer = -1 

desk.all_species = {}

print "\n--- Selected species: H0 "

for species in desk.dic_sel.keys():
    # species = random.choice(desk.dic_sel.keys())
    mat = desk.dic_sel[species]
    N = mat[1]
    selSpecies[species] = [N]
    
    ''' get only 4 elements 
    selSpecies[species].append(random.sample(list(range(N)), desk.sampleElems))
    
    get all:'''
    selSpecies[species].append([])
    selSpecies[species].append("H0")
    desk.all_species[species] = selSpecies[species]
          
    print species, selSpecies[species]
            
''' continue the counts to aggregate H0 species .... next Ha species '''

Lmincut = mat[5]
Lmaxmer = mat[6]

if opt == "mincut":
    desk.LH0 = Lmincut
else:
    desk.LH0 = Lmaxmer
    
    
print "Lmincut: %i Lmaxmer: %i bp"%(Lmincut, Lmaxmer)    

print "\n--- Randomizing Ha..."
notSelSpecies = {}
'''
# species = random.choice(desk.dic_not_sel.keys())
nSamp = int(len(desk.dic_not_sel.keys())*.8)
listSpecies = random.sample(desk.dic_not_sel.keys(), nSamp)
'''
curated = ["barutani", "virilis","mauritiana", "neocordata","equinoxialis","daruma",\
           "sechellia","kuntzei","emarginata"]
not_curated = ["lacertosa", "lutescens", "polychaeta","tsigana","prosaltans","biauraria","curvispina", "bipectinata",\
               "jambulina","lini","phalerata","rufa","affinis","eugracilis","auraria","serrata","quadraria","birchii",\
               "medioconstricta","funebris","transversa","nebulosa"]



print "\n--- Not selected species: Ha "
for species in desk.dic_not_sel.keys():
    if species in not_curated: continue
    
    mat = desk.dic_not_sel[species]
    N = mat[1]
    notSelSpecies[species] = [N]

    notSelSpecies[species].append([random.randint(0, N-1) for _ in range(5)])
    notSelSpecies[species].append("Ha")
    desk.all_species[species] = notSelSpecies[species]
    
    print species, notSelSpecies[species]

print "\nlooping"

desk.sample=False

if not desk.looping_data(selSpecies, hypothesis="H0", shuffle=shuffle):
    exit()
    
if not desk.looping_data(notSelSpecies, hypothesis="Ha"):
    exit()


desk.de_novo = False
desk.showGraph = True
desk.saveGraph = False
desk.saveData = False

desk.oneSample = False

desk.looping_JSD_all_samples()

        