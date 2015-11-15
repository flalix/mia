#-*- coding: utf-8 -*-
'''
Created on Sep 22, 2015

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica
'''
import os, platform
import random
import classes.Mutual_Information_Horizontal as MI
import numpy as np
import classes.Drosophila as dro
#import classes.Timing as crono
import classes.BioPythonClass as bioClass
import classes.JSD as JSD
#from scipy.stats.mstats import mquantiles
import classes.BarGraphic as graphPac

class Desktop():
    ''' we have:
            one sampled species: with **, e.g. pseudoobscura**
            nH0 species
            nHa species
            
            see: hyp_redo_hmi()
            
        here we 
        If analyzing the reference (one sample species) - maintain the sampled sequnce
        Otherise:
            Hypothesis can be H0 or Ha. 
            For H0:
                if species != reference (one sample species) - maintain the hole set
                    e.g. willistoni != pseudoobscura** 
                otherwise remove the sampled sequences
                    e.g. pseudoobscura == pseudoobscura**
            For Ha:
                if species != reference (one sample species) - maintain the hole set
                    e.g. neocordata != pseudoobscura** 
                otherwise remove the sampled sequences
                    e.g. mauritiana == mauritiana**
    '''     
    def __init__(self, filename_default, organism, gene, title, cutoffLength, cutoffNumSeq, showmessage=False):
        self.oneSample = False
        self.which = 'HMI'
        
        self.organism = organism
        self.seqType = "Gene"
        self.gene = gene 
        self.title = title
        self.set_gene_title()
        
        self.cutoffLength = cutoffLength
        self.cutoffNumSeq = cutoffNumSeq
        
        self.numOfLetters = 1
        self.withCorrection = True
        
        
        self.frame = 0
        self.offset = 0
        
        self.showmessage = showmessage

        self.platform = platform.system() # 'Linux'  Win32
        self.isWindows = (self.platform == 'Windows')

        self.home = os.path.expanduser("~")
        self.home = self.home.replace("\\","/")
        
        self.root = self.home + '/data/'          
        
        self.root_org   = self.home + '/data/%s/' %(self.organism)       
        self.root_files = self.home + '/data/%s/%s/' %(self.organism, self.gene_title)       
        self.rootFasta = self.home + '/data/%s/%s/fasta/'%(self.organism, self.gene_title)
        self.rootImage = self.home + '/data/%s/%s/image/'%(self.organism, self.gene_title)
        self.rootTable = self.home + '/data/%s/%s/files/'%(self.organism, self.gene_title)
        self.rootTree = self.home + '/data/%s/%s/trees/'%(self.organism, self.gene_title)
        self.rootEntropy = self.home + '/data/entropy/'

        rootTable = self.home + '/data/%s/%s/files/'%(organism, self.gene_title)
        self.filename_species_list =  rootTable + organism + '_' + self.gene_title + "_species_list.txt"  

        self.de_novo = False
        self.each_all = "each"
        self.mnat = True
        self.unit = 'mnat' if self.mnat else 'nat'
        self.factor = 1000 if self.mnat else 1
        self.norm = False
        
        self.numOfSEs = 2
        self.roundVal = 4
        
        self.refSpecies = None
        self.refSpeciesH0a = "H0"
       
        '''
        self.showGraph = False
        self.saveGraph = False
        self.saveData = True
        '''
        self.dpi = 300
        self.imgType = "png"
        self.isLog = False
        
        self.tk_root = None

        if self.organism == "Drosophila":
            self.dr = dro.Drosophila()
            
        self.entFile = bioClass.Entropy_File()
        self.gg = graphPac.GeneralGraphic(self) 


    def hyp_species_loop(self, nLoops=5):

        print "\n>>>H0 sampling"
        listH0_JSD = self.looping_H0_same_species(nLoops)

        print "\n>>>Ha sampling"
        listHa_JSD = self.looping_Ha_diff_species(nLoops)       

        return listH0_JSD, listHa_JSD


    def looping_H0_same_species(self, nLoops):
        lista = []
        self.jsd = JSD.Class_JSD(self)
        self.filename_correction = '_bias_corr' if self.withCorrection else ''
        
        self.hmi = MI.Mutual_Information_Horizontal(self, want_nuc_aa_list=False, want_shannon = False)
        if not self.hmi.ok:
            self.error_msg = 'Problems with Init_MI_Heat_Map_Hist_Graph.'
            return False
        
        
        for spec in sorted(self.all_species.keys()): # ["sturtevanti"]: # 
            '''  for each species simulate H0 '''
            # self.spRef = self.which_sp(self.organism, spec)
    
            ret, dList  = self.read_and_calc_H0(spec, nLoops=nLoops)
            if not ret:
                return None

            lista.append(dList)

        return lista
    
    
    def looping_Ha_diff_species(self, nLoops):
        lista = []
        self.jsd = JSD.Class_JSD(self)
        self.filename_correction = '_bias_corr' if self.withCorrection else ''
        self.entFile = bioClass.Entropy_File()
        
        self.hmi = MI.Mutual_Information_Horizontal(self, want_nuc_aa_list=False, want_shannon = False)
        if not self.hmi.ok:
            self.error_msg = 'Problems with Init_MI_Heat_Map_Hist_Graph.'
            return False
        
        specList = self.all_species.keys()
        n = len(specList)
        for i in range(n-1):
            for j in range(i+1, n):
                '''  for each pair of species simulate Ha '''

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

        self.hmi.ent.mySequence.set_seq_all_param(self, species, seqs=None,setSample=setSample)
        self.hmi.ent.mySequence.lenSequence = self.LH0

        # _ = self.time.start()
        k = .25 if N >= 20 else .5
        NSample = int(N * k)
        
        dList = []
        
        for jj in range(nLoops):
            stri = '%i/%i - %s'%(jj, nLoops, species)
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
                #print self.time.milli(), 'ms'
            
                if not self.withCorrection:           
                    arrayVal = np.array(ValList)
                    #arraySE = np.array(SElist)
                else:
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
        
        #list_distSE = []
        list_distMI = []
        
        setSample = False

        for jj in range(nLoops):
            stri = '%i/%i - %s x %s'%(jj, nLoops, species01, species02)
            self.desk.showmsg_obs(stri)           

            for species in [species01, species02]:
                N = self.all_species[species][0]
                knownSpecies = self.all_species[species][2]
                
                self.hmi.ent.mySequence.set_seq_all_param(self, species, seqs=None,setSample=setSample)
                self.hmi.ent.mySequence.lenSequence = self.LH0
    
                # _ = self.time.start()
                k = .25 if N >= 20 else .5
                NSample = int(N * k)
                
                dList = []
            
                sampledSeqs = random.sample(list(range(N)), NSample)

                ValList, ValCorrList, _, _, _ = \
                self.hyp_redo_hmi(species, True, N, sampledSeqs, knownSpecies)
                #_ = self.time.finish() 
                #print self.time.milli(), 'ms'

                if not self.withCorrection:           
                    arrayVal = np.array(ValList)
                    #arraySE = np.array(SElist)
                else:
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
            print "\n\n------------ Contingency  table -------------------\n"
            
            print "   |     True     |   False   |"
            print "---------------------------------"

            for key in self.contig_table.keys():
                m = self.contig_table[key]
                if key == "H0":
                    print "H0 |   TP = %i    |   FP = %i |"%(m[0], m[1])
                else:
                    print "Ha |   FN = %i    |   TN = %i |"%(m[0], m[1])
                print "---------------------------------"

    ''' must return: ValList, ValCorrList, SElist, SEcorrList '''
    def hyp_redo_hmi(self, species, isReference, N, sampledSeqs, knownSpecies):
         
        consensus = True if knownSpecies else False 
        numOfSeqs = len(sampledSeqs)
        
        hmi = MI.Mutual_Information_Horizontal(self, want_nuc_aa_list=False)
        
        if not hmi.ok:
            return None, None, None, None, None

        if knownSpecies:
            filename = ('%s_%s_%s_%s_%s_%iL_cutoff%i_consensus.fasta') %\
                 (self.organism, self.minmax, species, self.seqType, self.gene, self.cutoffLength, self.cutoffNumSeq)
        else:
            filename = ('%s_%s_%s_%s_%iL.fasta') %\
                 (self.organism, species, self.seqType, self.gene, self.cutoffLength)
            
        name = '%s %s %s'%(self.organism, species, self.gene_title)

        self.desk.showmsg_obs(">>reading %s"%(filename), same_line=True)
        
        if numOfSeqs == 0:
            self.error_msg = "numOfSeqs=%i, for %s"%(numOfSeqs,name)
            return None, None, None, None, None
        
        if not hmi.ent.mySequence.read_file_only_samples(self, filename, name, sampledSeqs, isReference):
            self.error_msg = "Could not find %s"%(filename)
            return None, None, None, None, None

        self.list_species_name_hmi = []
        
        seqs = []
        for rec in hmi.ent.mySequence.seqs:
            seqs.append(str(rec))

        ''' if sample = True we don't read and write files '''
        hmi.ent.mySequence.set_seq_all_param(self, species, seqs, setSample=False)
        
        sufix = hmi.ent.mySequence.sufix_hmi
        mat = [filename, name, species, hmi, sufix]
        self.list_species_name_hmi.append(mat)

        ''' always resample when Hypothesis Test '''
        ret, msg = hmi.prepareHorizontalMutualInfo(self, sufix, consensus=consensus, showmessage=self.showmessage, hyp_test=True)
            
        if not ret:
            self.showmsg_obs(msg)
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


    def set_gene_title(self):
        
        if self.gene != '':
            self.gene_title = self.gene
            self.graph_gene_title = self.gene
        else:
            self.gene_title = ''
            self.graph_gene_title = ''
        
        if self.title != '':
            if self.gene == '':
                self.gene_title = self.title
                self.graph_gene_title  = self.title
            else:
                self.gene_title += '_' + self.title
                self.graph_gene_title  += ' ' + self.title
                
        if self.gene_title == '':
            self.gene_title = 'undef'
            self.graph_gene_title = 'undef'
            
    def read_text(self, filename):
        try:
            if not os.path.isfile(filename):
                print 'Could not find %s'%(filename)
                return False, None
            
            f = open(filename, 'r')
            
            lines = f.readlines()
            for i in range(len(lines)):
                lines[i] = lines[i].rstrip('\n')
            
            fileError = False
        except:
            fileError = True
            
         
        if fileError:
            print 'Could not find %s'%(filename)
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


    def showmsg_obs(self, stri, same_line = False):
        if same_line: print stri,
        else: print stri

    
    def save_dic(self, filename, dics):
        self.entFile.write_file_data_dic2(filename, dics)
        
    def save_2lists(self, minmax, listH0, listHa):
        self.entFile.write_file(self.rootTable, "H0_%s.txt"%(minmax), str(listH0), showSaved=False)
        self.entFile.write_file(self.rootTable, "Ha_%s.txt"%(minmax), str(listHa), showSaved=False)

    def read_2lists(self):
        tudo_ok, listH0_JSD = self.entFile.read_file(self.rootTable, "H0_%s.txt"%(self.minmax))
        
        if tudo_ok:
            tudo_ok, listHa_JSD = self.entFile.read_file(self.rootTable, "Ha_%s.txt"%(self.minmax))
        else:
            listHa_JSD = []
            
        try:
            listH0_JSD = eval(listH0_JSD)
        except:
            print "could not convert list H0"
            listH0_JSD = []
            
        try:
            listHa_JSD = eval(listHa_JSD)
        except:
            print "could not convert list Ha"
            listHa_JSD = []    
            
            
        return tudo_ok, listH0_JSD,  listHa_JSD 
    
    
    def build_distributions_ROC(self, desk, listH0_JSD, listHa_JSD):
        
        if self.unit == "nat":
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
            

        # print len(listH0)
        # print len(listHa)
       
        self.contigency_table(self, listH0, listHa)
        
        self.gg.histogram_H0_Ha(self, listH0, listHa)
        
        pictureName = "ROC_JSD_HMI_%s"%(desk.minmax)
        self.gg.myPlot.print_graph(desk, self.gg.fig, pictureName, frame=desk.tk_root, stay=desk.stay)
        
        
    def contigency_table(self, desk, listH0, listHa):

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
            
            print "TP=%.3f, FP=%.3f, TN=%.3f, FN=%.3f"%(TP, FP, TN, FN)
            print "for %s, sensitivity=%.3f and specificity=%.3f"%(desk.minmax, sens, spec)
            print ""
            
            sensitivityList.append(sens)
            specificityList.append(spec)
 
        self.gg.ROC_Curve(desk, sensitivityList, specificityList)

        

