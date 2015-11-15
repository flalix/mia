'''
Created on Sep 14, 2015

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica

for each species calc if it is p-value < .05 on its JSD distribution
'''
import classes.BarGraphic as graphPac
gg = graphPac.GeneralGraphic()

import tests.Desktop_test as Desktop_test
import random
import classes.Drosophila as dro
 
organism = "Drosophila"      
gene = "Adh"      
title = ''
cutoffLength = 100
cutoffNumSeq = 10

filename_default = 'default.ini'
        
if organism == "Drosophila":
    dr = dro.Drosophila()
else:
    dr = None
    
desk = Desktop_test.Desktop(filename_default, organism, gene, title, cutoffLength, cutoffNumSeq)

opt = "mincut"
desk.minmax =  opt

desk.NminSamples = 6
desk.sampleElems = 4
#----- fix here ------------
desk.simulation = False
#---------------------------
shuffle = False
if shuffle: simulattion = False

print desk.filename_species_list

desk.build_listSpecies()
# print desk.cutoffNumSeq, desk.totList

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

desk.Lmincut = -1
desk.Lmaxmer = -1 

desk.all_species = {}
print "\n--- Selected species: H0 = %i species"%(len(desk.dic_sel.keys()))

for species in sorted(desk.dic_sel.keys()):
    # species = random.choice(desk.dic_sel.keys())
    mat = desk.dic_sel[species]
    N = mat[1]
    selSpecies[species] = [N]
    
    factor = .3 if N >= 20 else .5
    NSample = int(N * factor)

    selSpecies[species].append(random.sample(list(range(N)), NSample))
    selSpecies[species].append("H0")
    desk.all_species[species] = selSpecies[species]
          
    print species, NSample, selSpecies[species]


''' continue the counts to aggregate H0 species .... next Ha species '''

desk.Lmincut = mat[5]
desk.Lmaxmer = mat[6]

if opt == "mincut":
    desk.LH0 = desk.Lmincut
else:
    desk.LH0 = desk.Lmaxmer
    
    
print "Lmincut: %i Lmaxmer: %i bp"%(desk.Lmincut, desk.Lmaxmer)    

print "\n--- Randomizing Ha..."
notSelSpecies = {}
'''
# species = random.choice(desk.dic_not_sel.keys())
nSamp = int(len(desk.dic_not_sel.keys())*.8)
listSpecies = random.sample(desk.dic_not_sel.keys(), nSamp)

curated = ["barutani", "virilis","mauritiana", "neocordata","equinoxialis","daruma",\
           "sechellia","kuntzei","emarginata"]
not_curated = ["lacertosa", "lutescens", "polychaeta","tsigana","prosaltans","biauraria","curvispina", "bipectinata",\
               "jambulina","lini","phalerata","rufa","affinis","eugracilis","auraria","serrata","quadraria","birchii",\
               "medioconstricta","funebris","transversa","nebulosa"]
'''

    
factor = .5
    

for species in sorted(desk.dic_not_sel.keys()):
    # if species in not_curated: continue
    
    mat = desk.dic_not_sel[species]
    N = mat[1]
    notSelSpecies[species] = [N]
    
    NSample = int(N * factor)

    # notSelSpecies[species].append([random.randint(0, N-1) for _ in range(5)])
    notSelSpecies[species].append(random.sample(list(range(N)), NSample))
    notSelSpecies[species].append("Ha")
    desk.all_species[species] = notSelSpecies[species]
    
    print species, NSample, notSelSpecies[species]


print "\n--- Not selected curated species: Ha  = %i species"%(len(notSelSpecies.keys()))

create_samples = False

if create_samples:
    print "\nlooping: create samples / graphics / summary tables (selected and not sel)"
    
    desk.de_novo = False
    desk.showGraph = True
    desk.saveGraph = False
    desk.saveData = True
    
    desk.sample=False
    
    if not desk.looping_data(selSpecies, all_or_random = "all", hypothesis="H0",cmp_file="_sample_selected"):
        exit()
    
    if not desk.looping_data(notSelSpecies, all_or_random = "all", hypothesis="Ha",cmp_file="_sample_not_selected"):
        exit()
    
    print "\n---- end create_samples / looping selected and not selected ----"
    exit()
    
desk.de_novo = False
desk.showGraph = False
desk.saveGraph = False
desk.saveData = True


dics = {}
dics["mincut"] = []
dics["maxmer"] = []


lista = ["mincut","maxmer"]   
lista = ["mincut"]  
numLoops = 1

desk.oneSample = True
for minmax in lista:
    for i in range(numLoops):

        print ">>>>",i, minmax
        dic = desk.simulation_loop(minmax, showmessage=True)
        dics[minmax].append(dic)
        print ""

if desk.saveData:
    desk.save_dic(desk.rootTable+"summary_sens_spec.txt", dics)

sensitivityList, specificityList = [],[]

for minmax in lista:
    
    for dic in dics[minmax]:
        TP = dic["H0"][0]
        FP = dic["H0"][1]
        TN = dic["Ha"][1]
        FN = dic["Ha"][0]
        sens = TP / float(TP + FN)
        spec = TN / float(TN + FP)
        
        print "TP=%.3f, FP=%.3f, TN=%.3f, FN=%.3f"%(TP, FP, TN, FN)
        print "for %s, sensitivity=%.3f and specificity=%.3f"%(minmax, sens, spec)
        print ""
        
        sensitivityList.append(sens)
        specificityList.append(spec)
    

    gg.ROC_Curve(sensitivityList, specificityList)


'''
dics = []
dics.append(["mincut", {'H0': [11, 0], 'Ha': [1, 8]}])
dics.append(["maxmer", {'H0': [11, 0], 'Ha': [2, 7]}])

receiver operating characteristic (ROC)
is a graphical plot that illustrates the performance of a binary classifier  
system as its discrimination threshold is varied.

ROC:  TPR (y) x FPR (x)
      sensitivity x (1-specificity)
      
mincut  -  {'H0': [10, 1], 'Ha': [3, 6]}
TP=10.000, FP=1.000, TN=6.000, FN=3.000
for mincut, sensitivity=0.769 and specificity=0.857

mincut  -  {'H0': [10, 1], 'Ha': [2, 7]}
TP=10.000, FP=1.000, TN=7.000, FN=2.000
for mincut, sensitivity=0.833 and specificity=0.875

mincut  -  {'H0': [10, 1], 'Ha': [4, 5]}
TP=10.000, FP=1.000, TN=5.000, FN=4.000
for mincut, sensitivity=0.714 and specificity=0.833

mincut  -  {'H0': [10, 1], 'Ha': [2, 7]}
TP=10.000, FP=1.000, TN=7.000, FN=2.000
for mincut, sensitivity=0.833 and specificity=0.875

mincut  -  {'H0': [10, 1], 'Ha': [3, 6]}
TP=10.000, FP=1.000, TN=6.000, FN=3.000
for mincut, sensitivity=0.769 and specificity=0.857

maxmer  -  {'H0': [10, 1], 'Ha': [3, 6]}
TP=10.000, FP=1.000, TN=6.000, FN=3.000
for maxmer, sensitivity=0.769 and specificity=0.857

maxmer  -  {'H0': [10, 1], 'Ha': [2, 7]}
TP=10.000, FP=1.000, TN=7.000, FN=2.000
for maxmer, sensitivity=0.833 and specificity=0.875

maxmer  -  {'H0': [10, 1], 'Ha': [1, 8]}
TP=10.000, FP=1.000, TN=8.000, FN=1.000
for maxmer, sensitivity=0.909 and specificity=0.889

maxmer  -  {'H0': [10, 1], 'Ha': [2, 7]}
TP=10.000, FP=1.000, TN=7.000, FN=2.000
for maxmer, sensitivity=0.833 and specificity=0.875

maxmer  -  {'H0': [10, 1], 'Ha': [3, 6]}
TP=10.000, FP=1.000, TN=6.000, FN=3.000
for maxmer, sensitivity=0.769 and specificity=0.857 

x =[]
y=[]
sensitivity=0.769; specificity=0.857
x.append(sensitivity)
y.append(specificity)

sensitivity=0.833; specificity=0.875
x.append(sensitivity)
y.append(specificity)

sensitivity=0.714;  specificity=0.833
x.append(sensitivity)
y.append(specificity)

sensitivity=0.833;  specificity=0.875
x.append(sensitivity)
y.append(specificity)

sensitivity=0.769;  specificity=0.857
x.append(sensitivity)
y.append(specificity)

gg.ROC_Curve(x, y)

exit()    
'''

