'''
Created on Sep 22, 2015

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica

# H0: equal species - JSD --> 0
# Ha: different spcies - JSD >> 0

'''
import tests.Desktop_Hypothesis as Desktop_Hypothesis
 
organism = "Drosophila"      
gene = "Adh"      
title = ''
cutoffLength = 100
cutoffNumSeq = 10

filename_default = 'default.ini'
        

desk = Desktop_Hypothesis.Desktop(filename_default, organism, gene, title, cutoffLength, cutoffNumSeq)

desk.oneSample = False

desk.de_novo = False 
desk.showGraph = True
desk.saveGraph = True
desk.saveData = True
desk.stay = True

listaMinMax = ["mincut","maxmer"]   
listaMinMax = ["mincut"]  
listaMinMax = ["maxmer"]  


if not desk.de_novo:
    for minmax in listaMinMax:
        desk.minmax = minmax
        
        tudo_ok, listH0_JSD, listHa_JSD = desk.read_2lists()
        
        if tudo_ok:
            desk.build_distributions_ROC(desk, listH0_JSD, listHa_JSD)
            exit()
      
    
opt = listaMinMax[0]
desk.minmax =  opt

desk.NminSamples = 6
#----- fix here ------------
desk.simulation = False
shuffle = False
#---------------------------

print desk.filename_species_list

desk.build_listSpecies()
# print desk.cutoffNumSeq, desk.totList

desk.dic_sel = {}
desk.dic_not_sel = {}
desk.all_species = {}

for species in desk.dicParams.keys():
    mat = desk.dicParams[species]
    # print species, mat

    Nmincut = mat[2]
    Nmaxmer = mat[3]
    if opt == "mincut":
        N = Nmincut
    else:
        N = Nmaxmer

    '''
        [N, [], True]
        num of seqs, seqSels, known species
    '''
    if N >= desk.cutoffNumSeq:
        desk.dic_sel[species] = [N, [], True]
        desk.all_species[species] = [N, [], True]
    elif N >= desk.NminSamples:
        desk.dic_not_sel[species] = [N, [], False]
        desk.all_species[species] = [N, [], False]
        

desk.Lmincut = -1
desk.Lmaxmer = -1 

if opt == "mincut":
    desk.LH0 = desk.Lmincut
else:
    desk.LH0 = desk.Lmaxmer
    
desk.Lmincut = mat[5]
desk.Lmaxmer = mat[6]

print "Lmincut: %i Lmaxmer: %i bp"%(desk.Lmincut, desk.Lmaxmer)    

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



dics = {}
dics["mincut"] = []
dics["maxmer"] = []



for minmax in listaMinMax:
    print ">>>>",minmax
    listH0_JSD, listHa_JSD = desk.hyp_species_loop(nLoops=5)

    print ""

    if desk.saveData:
        desk.save_2lists(minmax, listH0_JSD, listHa_JSD)
    
    desk.build_distributions_ROC(minmax, listH0_JSD, listHa_JSD)

