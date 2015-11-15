#-*- coding: utf-8 -*-
'''
Created on 02/10/2015

@author: Flavio
'''
import os
import classes.Drosophila as dro
import random
import classes.Mutual_Information_Horizontal as MI
import numpy as np
import classes.Timing as crono
import classes.BioPythonClass as bioClass
import classes.JSD as JSD
from scipy.stats.mstats import mquantiles
import classes.BarGraphic as graphPac

class Classifier(object):
    '''
    Classifier: 
        build H0: same species
              Ha: different species
              calculates contingency tables
              calculates specificity and sensitivity
              calculares ROC
    '''


    def __init__(self, desk):
        self.desk = desk
        self.dr = desk.dr

        ''' ------------------------------------------------------- '''      
        self.entFile = bioClass.Entropy_File()
        self.gg = graphPac.GeneralGraphic(desk) 
        
        self.minmax = desk.minmax
        

    def hyp_species_loop(self, nLoops=5):

        self.desk.showmsg_obs("\n>>>Ha sampling")
        listHa_JSD = self.looping_Ha_diff_species(nLoops)       

        self.desk.showmsg_obs("\n>>>H0 sampling")
        listH0_JSD = self.looping_H0_same_species(nLoops)

        return listH0_JSD, listHa_JSD

    def looping_H0_same_species(self, nLoops):
        lista = []
        self.jsd = JSD.Class_JSD(self.desk)
        
        self.hmi = MI.Mutual_Information_Horizontal(self.desk, want_nuc_aa_list=False, want_shannon = False)
        if not self.hmi.ok:
            self.error_msg = 'Problems with Init_MI_Heat_Map_Hist_Graph.'
            return False
        
        
        for spec in sorted(self.all_species.keys()): # ["sturtevanti"]: # 
            '''  for each species simulate H0 '''
            # self.spRef = self.which_sp(self.desk.organism, spec)
    
            ret, dList  = self.read_and_calc_H0(spec, nLoops=nLoops)
            if not ret:
                return None

            lista.append(dList)

        return lista
    
    
    def looping_Ha_diff_species(self, nLoops):
        lista = []
        self.jsd = JSD.Class_JSD(self.desk)
        self.entFile = bioClass.Entropy_File()
        
        self.hmi = MI.Mutual_Information_Horizontal(self.desk, want_nuc_aa_list=False, want_shannon = False)
        if not self.hmi.ok:
            self.error_msg = 'Problems with Init_MI_Heat_Map_Hist_Graph.'
            return False
        
        specList = self.all_species.keys()
        n = len(specList)
        for i in range(n-1):
            for j in range(i+1, n):
                '''  for each pair of species simulate Ha '''

                if nLoops > 10:
                    nLoops = int(round(nLoops / 10.))

                ret, dList  = self.read_and_calc_Ha(specList[i], specList[j], nLoops=nLoops)
                if not ret:
                    return None
                
                lista.append(dList)

        return lista    


    ''' looping_JSD simple for H0 - species x itself '''
    def read_and_calc_H0(self,species, nLoops, showmessage=False):
        N = self.all_species[species][0]
        setSample = False

        ''' random.sample(list(range(N)), NSample) '''

        knownSpecies = self.all_species[species][2]

        self.hmi.ent.mySequence.set_seq_all_param(self.desk, species, seqs=None, which="HMI",setSample=setSample)
        self.hmi.ent.mySequence.lenSequence = self.LH0

        # _ = self.time.start()
        NSample = int(N * .5)
        
        dList = []
        
        for jj in range(nLoops):
            stri = '%i/%i - %s'%(jj+1, nLoops, species)
            self.desk.showmsg_obs(stri)            
              
            sampledSeqs = random.sample(list(range(N)), NSample)

            ''' stores a list of MI's for each species '''
            #list_distSE = []
            list_distMI = []
            
            for twice in [0,1]:
                isReference = True if twice == 0 else False

                ValList, ValCorrList, _, _, _ = \
                self.hyp_redo_hmi(species, isReference, N, sampledSeqs, knownSpecies)
    
                #_ = self.time.finish()
                #self.desk.showmsg_obs(self.time.milli(), 'ms')
            

                if not self.desk.withCorrection:
                    if not ValList:
                        continue
                             
                    arrayVal = np.array(ValList)
                    #arraySE = np.array(SElist)
                else:
                    if not ValCorrList:
                        continue
                             
                    arrayVal = np.array(ValCorrList)
                    #arraySE = np.array(SEcorrList)
    
                   
                ''' NORMALIZATION '''
                arrayVal /= np.sum(arrayVal)
                #arraySE *= arrayVal
                                                    
                list_distMI.append(arrayVal)
                #list_distSE.append(arraySE)
                ''' looping species + reference - saving distMI and distSE  ns.append(numOfSeqs)'''
        
            dList.append(self.jsd.calc_JS_Distance_simple(list_distMI[0], list_distMI[1]))

        return True, dList
    
    
    def read_and_calc_Ha(self,species01, species02, nLoops, showmessage=False):
        ''' stores a list of MI's for each pair of species '''
        ''' looping_JSD simple for Ha - species x another diff species '''
        

        setSample = False
        dList = []
        
        for jj in range(nLoops):
            stri = '%i/%i - %s x %s'%(jj, nLoops, species01, species02)
            self.desk.showmsg_obs(stri)
            
            #list_distSE = []
            list_distMI = []
            
            for species in [species01, species02]:
                N = self.all_species[species][0]
                knownSpecies = self.all_species[species][2]
                
                self.hmi.ent.mySequence.set_seq_all_param(self.desk, species, seqs=None, which="HMI", setSample=setSample)
                self.hmi.ent.mySequence.lenSequence = self.LH0

                NSample = int(N * .5)
           
                sampledSeqs = random.sample(list(range(N)), NSample)

                ValList, ValCorrList, _, _, _ = \
                self.hyp_redo_hmi(species, True, N, sampledSeqs, knownSpecies)
                #_ = self.time.finish() 
                #self.desk.showmsg_obs(self.time.milli(), 'ms')
            
                if not self.desk.withCorrection:           
                    if not ValList:
                        continue
                             
                    arrayVal = np.array(ValList)
                    #arraySE = np.array(SElist)
                else:
                    if not ValCorrList:
                        continue
                             
                    arrayVal = np.array(ValCorrList)
                    #arraySE = np.array(SEcorrList)
    
                   
                ''' NORMALIZATION '''
                arrayVal /= np.sum(arrayVal)
                #arraySE *= arrayVal
                           
                list_distMI.append(arrayVal)
                #list_distSE.append(arraySE)

                ''' looping 2 species '''
    
            dList.append(self.jsd.calc_JS_Distance_simple(list_distMI[0], list_distMI[1]))

        return True, dList

    def simulate_ROC(self, listH0_JSD, listHa_JSD, showmessage=False):
        if showmessage:
            self.desk.showmsg_obs("n------------ Contingency  table -------------------\n")
            
            self.desk.showmsg_obs("   |     True     |   False   |")
            self.desk.showmsg_obs("---------------------------------")

            for key in self.contig_table.keys():
                m = self.contig_table[key]
                if key == "H0":
                    self.desk.showmsg_obs("H0 |   TP = %i    |   FP = %i |"%(m[0], m[1]))
                else:
                    self.desk.showmsg_obs("Ha |   FN = %i    |   TN = %i |"%(m[0], m[1]))
                self.desk.showmsg_obs("---------------------------------")

    ''' must return: ValList, ValCorrList, SElist, SEcorrList '''
    def hyp_redo_hmi(self, species, isReference, N, sampledSeqs, knownSpecies):
    
        consensus = True if knownSpecies else False 
        numOfSeqs = len(sampledSeqs)
        
        hmi = MI.Mutual_Information_Horizontal(self.desk, want_nuc_aa_list=False)
        
        if not hmi.ok:
            return None, None, None, None, None

        if knownSpecies:
            filename = ('%s_%s_%s_%s_%s_%iL_cutoff%i_consensus.fasta') %\
                 (self.desk.organism, self.desk.minmax, species, self.desk.seqType, self.desk.gene_title, self.desk.cutoffLength, self.desk.cutoffNumSeq)
        else:
            filename = ('%s_%s_%s_%s_%iL.fasta') %\
                 (self.desk.organism, species, self.desk.seqType, self.desk.gene_title, self.desk.cutoffLength)
            
        name = '%s %s %s'%(self.desk.organism, species, self.desk.gene_title)

        if numOfSeqs == 0:
            self.desk.error_msg = "numOfSeqs=%i, for %s"%(numOfSeqs,name)
            return None, None, None, None, None

        if not hmi.ent.mySequence.read_file_only_samples(self.desk, filename, name, isReference, sampledSeqs, method=self.desk.rand_method_var.get()):
            self.desk.error_msg = "Could not find %s"%(filename)
            return None, None, None, None, None

        self.list_species_name_hmi = []
        
        seqs = []
        for rec in hmi.ent.mySequence.seqs:
            seqs.append(str(rec))

        ''' if sample = True we don't read and write files '''
        hmi.ent.mySequence.set_seq_all_param(self.desk, species, seqs, which="HMI", setSample=False)
        
        sufix = hmi.ent.mySequence.sufix_hmi
        mat = [filename, name, species, hmi, sufix]
        self.list_species_name_hmi.append(mat)

        ''' always resample when Hypothesis Test '''
        ret, msg = hmi.prepareHorizontalMutualInfo(self.desk, sufix, consensus=consensus, showmessage=self.desk.showmessage, hyp_test=True)
            
        if not ret:
            self.desk.showmsg_obs(msg)
            return None, None, None, None, None

        ''' end loop of selected species '''
        
        return hmi.arrayMI, hmi.arrayMIcorr, hmi.arraySE, hmi.arraySEcorr, hmi.arrayN


    def which_sp(self,organism,species,lim=4):
        if organism == 'Drosophila':
            sp = self.dr.mnemonic(species)
        else:
            lista = species.split("_")
            sp = ''
            for i in range(len(lista)):
                if i > 1:
                    lim = 2
                    
                if sp != '':
                    if i > 1:
                        sp += lista[i][:lim]
                    else:
                        sp += '_' + lista[i][:lim]
                else:
                    sp += lista[i][:lim]
        
        return sp



    def read_text(self, filename):
        try:
            if not os.path.isfile(filename):
                self.desk.showmsg_obs('Could not find %s'%(filename))
                return False, None
            
            f = open(filename, 'r')
            
            lines = f.readlines()
            for i in range(len(lines)):
                lines[i] = lines[i].rstrip('\n')
            
            fileError = False
        except:
            fileError = True
            
         
        if fileError:
            self.desk.showmsg_obs('Could not find %s'%(filename))
            return False, None
                
        return True, lines
    
    
    def read_text_to_species_list(self, filename):
        ret, lines = self.read_text(filename)
        
        if not ret:
            return None
        
        return lines
        
    def build_listSpecies(self):
        lista = self.read_text_to_species_list(self.filename_species_list)
        self.dicParams = {}
        self.totList = 0
    
        ''' two initial dummy lines with totals '''
        try:
            
            stri = lista[0]
            stri = stri.split(',')
            
            stri = lista[1]
            stri = stri.split(' ')
            #self.cutoffNumSeq = int(stri[1])
            
        except:
            pass
        
        finally:
            try:
                lista.pop(0)
                lista.pop(0)
            except:
                pass
        
        for line in lista:
            mat2 = line.split(",")
            species = mat2[1]
            ''' mat[0] = 'x' '''
            for i in range(2,8):
                mat2[i] = int(mat2[i])
                
            self.dicParams[species] = [mat2[0], mat2[2],mat2[3],mat2[4],mat2[5],mat2[6],mat2[7]]
            # mat = self.dicParams[species]
            self.totList += 1
            
        return self.cutoffNumSeq
        
     
    def find_species(self, desc):
        pos = desc.find("|")
        if pos < 1:
            return ""
        desc = desc[pos+1:]
        
        pos = desc.find("|")
        if pos < 1:
            return ""
        
        return desc[pos+1:].strip()

    
    def save_dic(self, filename, dics):
        self.entFile.write_file_data_dic2(filename, dics)
        
    def set_filenames(self):

        filenameH0 = "H0_%s%s%s.txt"%(self.minmax,self.desk.filename_correction,self.desk.stri_random)
        filenameHa = "Ha_%s%s%s.txt"%(self.minmax,self.desk.filename_correction,self.desk.stri_random)
        
        return filenameH0, filenameHa
        
    def save_2lists(self, listH0, listHa):
        filenameH0, filenameHa = self.set_filenames()

        
        self.entFile.write_file(self.desk.rootTable, filenameH0, str(listH0), showSaved=False)
        self.entFile.write_file(self.desk.rootTable, filenameHa, str(listHa), showSaved=False)

    def read_2lists(self):
        filenameH0, filenameHa = self.set_filenames()
        
        tudo_ok, listH0_JSD = self.entFile.read_file(self.desk.rootTable, filenameH0)
        
        if tudo_ok:
            tudo_ok, listHa_JSD = self.entFile.read_file(self.desk.rootTable, filenameHa)
        else:
            listHa_JSD = []
            
        try:
            listH0_JSD = eval(listH0_JSD)
        except:
            self.desk.showmsg_obs("could not convert list H0")
            listH0_JSD = []
            
        try:
            listHa_JSD = eval(listHa_JSD)
        except:
            self.desk.showmsg_obs("could not convert list Ha")
            listHa_JSD = []    
            
            
        return tudo_ok, listH0_JSD,  listHa_JSD 
    
    
    def build_distributions_ROC(self, listH0_JSD, listHa_JSD):
        
        if self.desk.unit == "nat":
            if len(listH0_JSD[0]) == 1:
                listH0 = [x[0] for x in listH0_JSD]
            else:
                listH0 = []
                for x in listH0_JSD:
                    for y in x:
                        listH0.append(y)
                
            if len(listH0_JSD[0]) == 1:
                listHa = [x[0] for x in listHa_JSD]
            else:
                listHa = []
                for x in listHa_JSD:
                    for y in x:
                        listHa.append(y)
        else:
            if len(listH0_JSD[0]) == 1:
                listH0 = [x[0]*1000 for x in listH0_JSD]
            else:
                listH0 = []
                for x in listH0_JSD:
                    for y in x:
                        listH0.append(y*1000)
                
            if len(listH0_JSD[0]) == 1:
                listHa = [x[0]*1000 for x in listHa_JSD]
            else:
                listHa = []
                for x in listHa_JSD:
                    for y in x:
                        listHa.append(y*1000)
            

        # self.desk.showmsg_obs(len(listH0))
        # self.desk.showmsg_obs(len(listHa))
       
        self.contigency_table(listH0, listHa)
        
        self.gg.histogram_H0_Ha(self.desk, listH0, listHa, self.desk.label_random)
        
     
        pictureName = "ROC_JSD_HMI_%s%s%s"%(self.minmax,self.desk.filename_correction,self.desk.stri_random)
        self.gg.myPlot.print_graph(self.desk, self.gg.fig, pictureName, frame=self.desk.tk_root, stay=True)
        
        
    def contigency_table(self, listH0, listHa):

        mean0 = np.mean(listH0)
        meanA = np.mean(listHa)
        
        sensitivityList, specificityList = [],[]

        
        mat = np.linspace(mean0, meanA, num=50)
        
        for i in range(len(mat)):
            cutoff = mat[i]
            TP = np.sum([1 for x in listH0 if x <= cutoff])
            FP = np.sum([1 for x in listH0 if x >  cutoff])
            TN = np.sum([1 for x in listHa if x >= cutoff])
            FN = np.sum([1 for x in listHa if x <  cutoff])
            
            sens = TP / float(TP + FN)
            spec = TN / float(TN + FP)
            
            '''
            self.desk.showmsg_obs("TP=%.3f, FP=%.3f, TN=%.3f, FN=%.3f"%(TP, FP, TN, FN))
            self.desk.showmsg_obs("for %s, sensitivity=%.3f and specificity=%.3f"%(desk.minmax, sens, spec))
            self.desk.showmsg_obs("")
            '''
            
            sensitivityList.append(sens)
            specificityList.append(spec)
 
        self.gg.ROC_Curve(self.desk, sensitivityList, specificityList, self.desk.label_random)



    '''-------------- end --------------------'''
   
   
    def looping_JSD_one_sample(self):
        
        self.jsd = JSD.Class_JSD(self.desk)

        ''' -----------------------------------------'''  
        
            
        self.time = crono.Timing()
        self.entFile = bioClass.Entropy_File()
        self.set_my_cluster_filenames()
        
        name = '%s %s %s'%(self.desk.organism, self.gene, self.title)
        self.names = [[self.cluster_input_filename, name]]

        self.hmi = MI.Mutual_Information_Horizontal(self.desk, want_nuc_aa_list=False, want_shannon = False)
        if not self.hmi.ok:
            self.error_msg = 'Problems with Init_MI_Heat_Map_Hist_Graph.'
            return False

        if self.desk.organism == "Drosophila":
            self.dr = dro.Drosophila()
        
        # filename_dist_mat = self.cluster_input_filename
                   
        ''' always read and calc in looping_JSD_one_sample'''     
        ret, dicMI, _, _  = self.read_and_calc_one_sample()
       
        if not ret:
            return False
        
        xlabel = 'species x species'
        
        ylabel = 'JSD(%s)'%(self.desk.unit)

        if self.desk.showGraph or self.desk.saveGraph:
            self.jsd.plot_JSD_MI(self.desk, dicMI, self.desk.filename_correction, xlabel, ylabel, onlyRef=True, stay=True)

        
        return True
    def set_my_cluster_filenames(self):
        self.desk.set_gene_title()
        
        if self.frame == 0:
            self.str_frame1 = ''
            self.str_frame2 = ''
        else:
            self.str_frame1 = 'frame=%i_'%(self.frame)
            self.str_frame2 = ' frame=%i, '%(self.frame)

        self.sufixTitle = '%s %s%s %s, gene=%s%s NOL=%i' %\
                (self.desk.organism, self.minmax, self.desk.str_correction, self.seqType, self.desk.graph_gene_title, self.str_frame2, self.numOfLetters)
        self.sufix =     ('%s_%s_%s_%s_%sNOL%i_%iL_cutoff%i') %\
              (self.desk.organism, self.minmax, self.seqType, self.gene_title, self.str_frame1, self.numOfLetters, self.cutoffLength, self.cutoffNumSeq)

        self.prefix = 'JSD_HMI_Sample_'
        

        self.title_jsd = 'JSD(Horizontal Mutual Information) %s - %s\n%s'%\
            (self.desk.unit, self.refSpecies, self.sufixTitle)


        self.cluster_input_filename = 'JSD_HMI_Sample_%s%s.txt' % (self.sufix, self.desk.filename_correction)
        self.striMatrix = "Horizontal MI distance matrix"

    
    ''' looping_JSD_one_sample calls '''
    def read_and_calc_one_sample(self, showmessage=False):
        seqSpecies = []
        ''' stores a list of MI's for each species '''
        list_distSE = []
        list_distMI = []
        ns = []
       
        listSpecies = []
        listSpecJSD = []
        
        speciesList = sorted(self.all_species.keys())
        ''' the 0 is the reference ... all others are all_species to be compared '''
        k = len(self.all_species.keys()) + 1
        
        for iLoop in range(k):
            if iLoop == 0:
                species = self.refSpecies
                setSample = True
                sampleStr = "**"
                
                isReference = True
                N = self.N
                sampledSeqs = self.sampledSeqs
                hypothesis = self.refSpeciesH0a
            else:
                species = speciesList[iLoop-1]
                setSample = False
                sampleStr = ""

                isReference = False
                mat = self.all_species[species]
                N = mat[0]
                sampledSeqs = mat[1]
                hypothesis = mat[2]
                
            if isReference:
                stri = '\n%i/%i) %s #Seq=%i selected %s'%(iLoop+1, len(self.all_species.keys())+1, species, N, str(sampledSeqs))
                numOfSeqs = len(sampledSeqs)
            elif species == self.refSpecies:
                stri = '%i/%i) %s #Seq=%i complement of %s'%(iLoop+1, len(self.all_species.keys())+1, species, N, str(sampledSeqs))
                numOfSeqs = N - len(sampledSeqs)
            else:
                stri = '%i/%i) %s #Seq=%i all sequences'%(iLoop+1, len(self.all_species.keys())+1, species, N)
                numOfSeqs = N

            self.desk.showmsg_obs(stri)
            
            sp = self.which_sp(self.desk.organism, species)

            listSpecJSD.append(sp+sampleStr)

            L = self.LH0

            listSpecies.append(sp+sampleStr)

            if numOfSeqs == 0 or L == 0:
                self.error_msg = "numOfSeqs=%i, L=%i, for %s"%(numOfSeqs, L, species)
                continue
                          
            seqSpecies.append([species,sp,numOfSeqs])
            ns.append(numOfSeqs)

            self.hmi.ent.mySequence.set_seq_all_param(self.desk, species, seqs=None, which="HMI", setSample=setSample)
            self.hmi.ent.mySequence.lenSequence = L
            
            sufix = self.hmi.ent.mySequence.sufix_hmi
            filenameMI = 'HMI_%s_mi.txt'%(sufix) # sampleStr
                            
            # self.desk.showmsg_obs('>> %s sample len = %i ***** ' % (filenameMI, numSeqs))
            ''' hmi_Drosophila_americana_Gene_frame0_NOL1_400L_cutoff10_mij.txt '''
            if not self.entFile.setNames(self.desk.rootTable, filename=filenameMI):
                self.showmsg_obs('Could not find %s'%(self.desk.rootTable))
                return False, None, None, None

            # _ = self.time.start()

            if iLoop == 0 or species == self.refSpecies:
                ValList, ValCorrList, SElist, SEcorrList, _ = self.redo_hmi(species, isReference, N, sampledSeqs, hypothesis)
            else:
                ok, _,  ValList, ValCorrList, _ = \
                self.entFile.read_HMI_files(root=self.desk.rootTable, sufix=sufix, which='mi' ,showKeys=self.showmessage)
    
                if not ok:
                    self.desk.showmsg_obs("recalc: redo_hmi")
                    ValList, ValCorrList, SElist, SEcorrList, _ = self.redo_hmi(species, isReference, N, sampledSeqs, hypothesis)
                else:
                    ok, _,  SElist, SEcorrList, _ = \
                    self.entFile.read_HMI_files(root=self.desk.rootTable, sufix=sufix, which='se', showKeys=self.showmessage)
        
                    if not ok:
                        self.desk.showmsg_obs("recalc: redo_hmi")
                        ValList, ValCorrList, SElist, SEcorrList, _ = self.redo_hmi(species, isReference, N, sampledSeqs, hypothesis)
           
                    
            #_ = self.time.finish()
            #self.desk.showmsg_obs(self.time.milli(), 'ms')
        
            if not self.desk.withCorrection:           
                arrayVal = np.array(ValList)
                arraySE = np.array(SElist)
            else:
                arrayVal = np.array(ValCorrList)
                arraySE = np.array(SEcorrList)

               
            ''' NORMALIZATION '''
            arrayVal /= np.sum(arrayVal)
            arraySE *= arrayVal
                                                
            list_distMI.append(arrayVal)
            list_distSE.append(arraySE)
            
            iLoop += 1
            ''' looping species + reference - saving distMI and distSE  ns.append(numOfSeqs)'''


        dicMI, striJSD, sError = self.jsd.calc_JSD(listSpecies, list_distMI, list_distSE, parSE=self.numOfSEs,unit=self.desk.unit)


        self.desk.result = {}
        vals = []
        d = -1
        strRef = ""
        for key in dicMI.keys():
            if key.find("**") > 0:
                v = dicMI[key][0]
                vals.append(v)
                
                elems = key.split(" x ")
                if elems[0] == self.spRef or elems[1] == self.spRef:
                    d = v
                    strRef = key
             
        
        mean = np.mean(vals)*self.desk.factor
        self.desk.result["mean"] = round(mean,1)
        
        pvalue_inf, median = mquantiles(vals, prob=[0.05, 0.5])*self.desk.factor
        self.desk.result["pvalue_inf"] = round(pvalue_inf,1)
        self.desk.result["median"] = round(median,1) 

        # self.desk.showmsg_obs("\n\nSelected sample distribution -------------------")
        # self.desk.showmsg_obs("mean=%.1f, median=%.1f, pvalue_inf=%.1f"%(self.desk.result["mean"], self.desk.result["median"], self.desk.result["pvalue_inf"]))
        # self.desk.showmsg_obs("------------------------------------------------")
        
        if d == -1:
            stri = "Nothing found."
        else:
            d = round(d*self.desk.factor,1)
            
            if self.refSpeciesH0a not in self.contig_table.keys():
                self.contig_table[self.refSpeciesH0a] = [0, 0]
                
            if d <= pvalue_inf:
                if self.refSpeciesH0a == "H0":
                    self.contig_table[self.refSpeciesH0a][0] += 1
                else:
                    self.contig_table[self.refSpeciesH0a][1] += 1
                stri = "%s is statistically significant with d=%.1f and p-value=%.1f"%(strRef, d, self.desk.result["pvalue_inf"])
            else:
                if self.refSpeciesH0a == "H0":
                    self.contig_table[self.refSpeciesH0a][1] += 1
                else:
                    self.contig_table[self.refSpeciesH0a][0] += 1
                stri = "%s is not statistically significant with d=%.1f and p-value=%.1f"%(strRef, d, self.desk.result["pvalue_inf"])

        self.desk.showmsg_obs(stri)
        self.desk.result["result"] = stri
        
        return True, dicMI, striJSD, sError

        
    ''' must return: ValList, ValCorrList, SElist, SEcorrList '''
    def redo_hmi(self, species, isReference, N, sampledSeqs, hypothesis):
         
        consensus = True if hypothesis == "H0" else False 
        numOfSeqs = len(sampledSeqs)
        
        hmi = MI.Mutual_Information_Horizontal(self.desk, want_nuc_aa_list=False)
        
        if not hmi.ok:
            return None, None, None, None, None

        if hypothesis == "Ha":
            filename = ('%s_%s_%s_%s_%iL.fasta') %\
                 (self.desk.organism, species, self.seqType, self.gene, self.cutoffLength)
        else:
            filename = ('%s_%s_%s_%s_%s_%iL_cutoff%i_consensus.fasta') %\
                 (self.desk.organism, self.minmax, species, self.seqType, self.gene, self.cutoffLength, self.cutoffNumSeq)

        
        name = '%s %s %s'%(self.desk.organism, species, self.gene_title)

        if numOfSeqs == 0:
            self.error_msg = "numOfSeqs=%i, for %s"%(numOfSeqs,name)
            return None, None, None, None, None
        
        if not hmi.ent.mySequence.read_file_one_sample(self.desk, filename, name, isReference, species, N, self.sampledSeqs, hypothesis):
            self.error_msg = "Could not find %s"%(filename)
            return None, None, None, None, None

        self.list_species_name_hmi = []
        
        seqs = []
        for rec in hmi.ent.mySequence.seqs:
            seqs.append(str(rec))

        ''' if sample = True we don't read and write files '''
        hmi.ent.mySequence.set_seq_all_param(self.desk, species, seqs, which="HMI", setSample=False)
        
        sufix = hmi.ent.mySequence.sufix_hmi
        mat = [filename, name, species, hmi, sufix]
        self.list_species_name_hmi.append(mat)

        ''' resample Ha, get rid from sampledSeqs in H0 '''
        if species == self.refSpecies:
            ret, msg = hmi.prepareHorizontalMutualInfo_OneSample(self.desk, sampledSeqs, sufix=sufix, consensus=consensus, LH0a=self.LH0, showmessage=self.showmessage)
        else:
            ret, msg = hmi.prepareHorizontalMutualInfo(self.desk, sufix, consensus=consensus, LH0a=self.LH0, showmessage=self.showmessage)
            
        if not ret:
            self.showmsg_obs(msg)
            return None, None, None, None, None

        ''' end loop of selected species '''
        
        return hmi.arrayMI, hmi.arrayMIcorr, hmi.arraySE, hmi.arraySEcorr, hmi.arrayN

    def open_distance_matrix_file(self, filename):
            #open the file assuming the data above is in a file called 'dataFile'
            if not os.path.exists(filename):
                return False, None, None, None
                
            # self.desk.showmsg_obs('opening %s'%(filename))
        
            try:
                inFile = open(filename,'r')
                #save the column/row headers (conditions/genes) into an array
                colHeaders = inFile.next().strip().split()[1:]
                rowHeaders = []
                dataMatrix = []
                
                for line in inFile:
                    data = line.strip().split('\t')
                    rowHeaders.append(data[0])
                    dataMatrix.append([float(x) for x in data[1:]])
            except:
                return False, None, None, None
        
            return True, rowHeaders, colHeaders, dataMatrix
        
    def simulation_loop(self, minmax, showmessage=False):
        self.minmax =  minmax

        if minmax == "mincut":
            self.LH0 = self.Lmincut
        else:
            self.LH0 = self.Lmaxmer        

        self.contig_table = {}

        i = 0
        for spec in sorted(self.all_species.keys()):  # ["triauraria"]:  
            ''' for each species calc if it is p-value < .05 on its JSD distribution '''
            selSpec = self.all_species[spec]
            self.refSpecies = spec
            self.desk.refSpecies = spec
            self.N = selSpec[0]
            self.sampledSeqs = selSpec[1]
            self.refSpeciesH0a = selSpec[2]
            
            self.spRef = self.which_sp(self.desk.organism, spec)

            # self.showGraph = True if i <= 3 else False

            self.looping_JSD_one_sample()
            i += 1

        
                   
        if showmessage:
            self.desk.showmsg_obs("\n\n------------ Contingency  table -------------------\n")
            
            self.desk.showmsg_obs("   |     True     |   False   |")
            self.desk.showmsg_obs("---------------------------------")

            for key in self.contig_table.keys():
                m = self.contig_table[key]
                if key == "H0":
                    self.desk.showmsg_obs("H0 |   TP = %i    |   FP = %i |"%(m[0], m[1]))
                else:
                    self.desk.showmsg_obs("Ha |   FN = %i    |   TN = %i |"%(m[0], m[1]))
                self.desk.showmsg_obs("---------------------------------")
                
        return self.contig_table
