'''
Created on Sep 14, 2015
Update  on Sep 17, 2015 - redo_hmi() + docs

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica
'''
import os, platform
import matplotlib.pyplot as plt
import gc
import classes.Mutual_Information_Horizontal as MI
import classes.BarGraphic as graphPac
import numpy as np
import classes.Drosophila as dro
import classes.Timing as crono
import classes.BioPythonClass as bioClass
import classes.JSD as JSD
from scipy.stats.mstats import mquantiles

class Desktop():
    ''' we have:
            one sampled species: with **, e.g. pseudoobscura**
            nH0 species
            nHa species
            
            see: redo_hmi()
            
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
                
          
    def looping_data(self, selSpecies, all_or_random="random", hypothesis="H0", shuffle=False, cmp_file=""):
        print("\n---------------------------------")
        plt.close("all")
        plt.clf()
                    
        self.list_species_name_hmi = []
        
        iLoop = 0
        
        for species in selSpecies.keys():
            mat = selSpecies[species]
            N = mat[0]
            sampledSeqs = mat[1]
            H0a = mat[2]
            consensus = True if H0a == "H0" else False
            setSample = True if sampledSeqs else False
            print "From %s selected %s in %i seqs."%(species, str(sampledSeqs), N)
            
            hmi = MI.Mutual_Information_Horizontal(self, want_nuc_aa_list=False)
            
            if not hmi.ok:
                return False
            
            iLoop += 1
    
            numOfSeqs = len(sampledSeqs)
                
            stri = '%i/%i) %s #Seq=%i'%(iLoop, len(selSpecies.keys()), species, numOfSeqs)
    
            if hypothesis and hypothesis == "Ha":
                filename = ('%s_%s_%s_%s_%iL.fasta') %\
                     (self.organism, species, self.seqType, self.gene, self.cutoffLength)
            else:
                filename = ('%s_%s_%s_%s_%s_%iL_cutoff%i_consensus.fasta') %\
                     (self.organism, self.minmax, species, self.seqType, self.gene, self.cutoffLength, self.cutoffNumSeq)

            print stri, filename
            name = '%s_%s%s'%(self.organism, species, self.gene_title)

            if numOfSeqs == 0:
                self.error_msg = "numOfSeqs=%i, for %s"%(numOfSeqs,name)
                continue
            
            if self.simulation:
                hmi.ent.mySequence.fasta_file_simulation(self, L=self.LH0, N=self.sampleElems)
                '''
                print hmi.ent.mySequence.numSequences
                print hmi.ent.mySequence.lenSequence
                for k in range(self.Nsamples):
                    print hmi.ent.mySequence.seqs[k]
                '''
            else:
                if not hmi.ent.mySequence.read_file_new(self, filename, name, False, all_or_random, sampledSeqs, shuffle):
                    self.error_msg = "Could not find %s"%(filename)
                    return False
    
       
            seqs = []
            for rec in hmi.ent.mySequence.seqs:
                seqs.append(str(rec))

            
            '''  if simulation sim_ + species'''
            str_sim = "sim_" if self.simulation else ""
            hmi.ent.mySequence.set_seq_all_param(self, str_sim+species, seqs, setSample=setSample)
            
            sufix = hmi.ent.mySequence.sufix_hmi
            mat = [filename, name, species, hmi, sufix]
            self.list_species_name_hmi.append(mat)
            
            ret, msg = hmi.prepareHorizontalMutualInfo(self, sufix=sufix, consensus=consensus, LH0a=self.LH0, showmessage=self.showmessage)

            if not ret:
                self.showmsg_obs(msg)
                return False

        ''' end loop of selected species '''
    
        self.save_tables_show_graph(cmp_file) 
        ''' end of looping_data '''
        return True 

    ''' must return: ValList, ValCorrList, SElist, SEcorrList '''
    def redo_hmi(self, species, isReference, N, sampledSeqs, hypothesis):
         
        consensus = True if hypothesis == "H0" else False 
        numOfSeqs = len(sampledSeqs)
        
        hmi = MI.Mutual_Information_Horizontal(self, want_nuc_aa_list=False)
        
        if not hmi.ok:
            return None, None, None, None, None

        if hypothesis == "Ha":
            filename = ('%s_%s_%s_%s_%iL.fasta') %\
                 (self.organism, species, self.seqType, self.gene, self.cutoffLength)
        else:
            filename = ('%s_%s_%s_%s_%s_%iL_cutoff%i_consensus.fasta') %\
                 (self.organism, self.minmax, species, self.seqType, self.gene, self.cutoffLength, self.cutoffNumSeq)

        
        name = '%s %s %s'%(self.organism, species, self.gene_title)


        if numOfSeqs == 0:
            self.error_msg = "numOfSeqs=%i, for %s"%(numOfSeqs,name)
            return None, None, None, None, None
        
        if not hmi.ent.mySequence.read_file_one_sample(self, filename, name, isReference, species, N, self.sampledSeqs, hypothesis):
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

        ''' resample Ha, get rid from sampledSeqs in H0 '''
        if species == self.refSpecies and self.oneSample:
            ret, msg = hmi.prepareHorizontalMutualInfo_OneSample(self, sampledSeqs, sufix=sufix, consensus=consensus, LH0a=self.LH0, showmessage=self.showmessage)
        else:
            ret, msg = hmi.prepareHorizontalMutualInfo(self, sufix, consensus=consensus, LH0a=self.LH0, showmessage=self.showmessage)
            
        if not ret:
            self.showmsg_obs(msg)
            return None, None, None, None, None

        ''' end loop of selected species '''
    
        # self.save_tables_show_graph()
        
        return hmi.arrayMI, hmi.arrayMIcorr, hmi.arraySE, hmi.arraySEcorr, hmi.arrayN


    def save_tables_show_graph(self,cmp_file=""):
        bias_corr_list = [self.withCorrection]

        print("-----------------------------------------------------------")
        for bias_corr in bias_corr_list:
            iLoop = 0
            speciesParams = []
            listMI_Anova = []       
                 
            for mat in self.list_species_name_hmi:
                iLoop += 1
                filename, name, species, hmi, sufix = mat
                

                numOfSeqs = self.sampleElems
                L = self.LH0

   
                print("%i) %s - %s, %s"%(iLoop, name, self.minmax, "bias corr." if bias_corr else 'no corr.' ))
    
                sp = self.which_sp(self.organism, species)
                        
                self.withCorrection = bias_corr
                
                if self.withCorrection:
                    str_correction = ' (bias corr.)'
                    filename_correction = '_bias_corr'
                else:
                    str_correction = ''
                    filename_correction = ''        

                if self.frame == 0:
                    sFrame = ''
                else:
                    sFrame = ', frame=' + str(self.frame)
 
                title = 'HMI: %s %s, %s %s\n%s%s; %i seqs; len=%i%s; #letter=%i'%\
                        (self.organism, species, self.seqType, self.gene, self.minmax, str_correction, numOfSeqs, L, sFrame, self.numOfLetters)

                                
                ''' HMI mean and sdv for each k from n sequences '''
                if not self.withCorrection:           
                    arrayMIthis = np.array(hmi.arrayMI)
                    arraySE     = np.array(hmi.arraySE)
                else:
                    arrayMIthis = np.array(hmi.arrayMIcorr)
                    arraySE     = np.array(hmi.arraySEcorr)
    
               
                if self.mnat:
                    arrayMIthis *= 1000
                    arraySE *= 1000
    
                    roundVal = 1
                    self.unit = 'mnat'
                else:
                    roundVal = 4
                    self.unit = 'nat'
                    
                if self.norm:
                    arrayMIthis /= self.numOfLetters
                    arraySE /= self.numOfLetters
                
                # self.showmsg_obs('  >>> species %s with %i sequences'%(species, numOfSeqs))
                '''  first time build abcissa
                         if choose == 1 and self.frame == 1
                         or self.frame == 0
                '''
                xSeq = []
                
                if self.frame < 2:
                    for k in range( len(arrayMIthis)):
                        xSeq.append(k+3)
                else:
                    ''' frame1: k = {3,6,9, ... 3n} '''
                    for k in range( len(arrayMIthis)):
                        xSeq.append((k+1)*3)
                                        
           
                arrayVal = []
                ySeqSEFinal = []


                ''' if only one curve '''
                if self.frame < 2:
                    arrayVal = arrayMIthis
                    ySeqSEFinal = arraySE
                else:
                    ''' if choose > 1 there are 3 others curves, and k jumps 3 x 3 '''
                    for m in range(len(arrayMIthis)):
                        arrayVal.append(arrayMIthis[m])
                        arrayVal.append(arrayMIthis[m])
                        arrayVal.append(arrayMIthis[m])
                        
                        ySeqSEFinal.append(arraySE[m])
                        ySeqSEFinal.append(arraySE[m])
                        ySeqSEFinal.append(arraySE[m])
    
                self.meanY = np.round(np.mean(arrayVal), roundVal)
                self.medianY = np.round(np.median(arrayVal), roundVal)
                self.stdY  = np.round(np.std(arrayVal), roundVal)
                self.maxY = np.round(np.max(arrayVal), roundVal)
                self.minY = np.round(np.min(arrayVal), roundVal)
                
                            
                speciesParams.append([species,sp,numOfSeqs, self.meanY, self.medianY, self.stdY, self.maxY, self.minY])
                listMI_Anova.append(arrayVal)
                
                if self.showGraph or self.saveGraph:
                    if self.frame < 2:
                        gHist = graphPac.Histogram_FreqDistribution(self, title)
                        
                        gHist.meanY = self.meanY
                        gHist.medianY = self.medianY
                        gHist.stdY  = self.stdY
                        gHist.maxY = self.maxY
                        gHist.minY = self.minY
                        gHist.desk = self
                                        
                        gHist.plot_H_MI('V', xSeq=xSeq, ySeq=arrayVal, ySE=arraySE, showError=True, \
                                        unit=self.unit, roundVal=roundVal)                
                        
                        gHist.plot_H_MI('H', xSeq=xSeq, ySeq=arrayVal, ySE=ySeqSEFinal, showError=True, \
                                        unit=self.unit, roundVal=roundVal)
                        
                        gHist.densityBar("H", seq=arrayVal, unit=self.unit, roundVal = roundVal)
                    else:
                        gHist.sameBar(xSeq=xSeq, arrayMIthis=arrayVal,linestyleCode=self.arrLinestyleCode[self.frame], color=self.arrColor[self.frame])
                            
                    
                    if (self.frame == 0) or (self.frame==3): 
                        ''' sample in filename '''   
                        pictureName = 'HMI_%s%s%s'%(sufix, filename_correction, cmp_file)
                        gHist.myPlot.print_graph(self, gHist.fig, pictureName, frame=self.tk_root, stay=True)
    
                    '''
                    del gHist
                    plt.cla()
                    plt.clf()
                    gc.collect()
                    '''

            if self.saveData :
                sufix = ('%s_%s_%s_%s_frame%i_NOL%i_%iL_cutoff%i') %\
                        (self.organism, self.minmax, self.seqType, self.gene, self.frame, self.numOfLetters, self.cutoffLength, self.cutoffNumSeq)

                ''' summary with sample in filename '''
                filename = 'HMI_summary_%s%s%s.txt' % (sufix, filename_correction, cmp_file)        
                stri = hmi.ent.calc_anova(listMI_Anova, sufix)
        
                hmi.ent.print_data_summary(self, speciesParams, roundVal=roundVal, filename=filename, stri=stri)

        self.failure = False
        self.error_msg = 'Task ended. All right.'
        
        ''' in the end of save_tables_show_graph clean memory '''
        del self.list_species_name_hmi
        
        try:
            gc.collect() 
        except:
            pass  
                
        return True

    def looping_JSD_all_samples(self):
        plt.close("all")
        plt.clf()

        self.jsd = JSD.Class_JSD(self)

        self.filename_correction = '_bias_corr' if self.withCorrection else ''
                      
        if self.withCorrection:
            self.filename_correction = '_bias_corr'
        else:
            self.filename_correction = ''
            
        self.time = crono.Timing()
        
        name = '%s %s %s'%(self.organism, self.gene, self.title)

        self.entFile = bioClass.Entropy_File()

        self.set_my_cluster_filenames()
        
        self.names = [[self.cluster_input_filename, name]]

        self.hmi = MI.Mutual_Information_Horizontal(self, want_nuc_aa_list=False, want_shannon = False)
        if not self.hmi.ok:
            self.error_msg = 'Problems with Init_MI_Heat_Map_Hist_Graph.'
            return False

        if self.organism == "Drosophila":
            self.dr = dro.Drosophila()
        
        filename_dist_mat = self.cluster_input_filename
                        
        if not os.path.exists(self.rootTable +  filename_dist_mat) or self.de_novo:
            ''' if not found in looping_JSD_all_samples read and calc '''
            ret, dicMI, stri, sError  = self.read_and_calc_all_samples()
        else:
            ret, dicMI, stri, sError = self.read_distance_matrix(filename_dist_mat)
        
        if not ret:
            return False
        
        xlabel = 'species x species'
        
        ylabel = 'JSD(%s)'%(self.unit)

        if self.saveData:
            self.showmsg_obs('\n ============ %s ============'%(self.striMatrix))
            self.showmsg_obs(stri)
            self.showmsg_obs('============ SE %s ============'%(self.striMatrix))
            self.showmsg_obs(sError)
            self.showmsg_obs('============ %s ============\n'%(self.striMatrix))

            if not self.hmi.entFile.write_file(self.rootTable, filename_dist_mat, stri, True):
                return False
            
            filename = filename_dist_mat.replace('.txt','') + '_se.txt'
            if not self.hmi.entFile.write_file(self.rootTable, filename, sError):
                return False

            striSummary = self.hmi.print_MI_species_data(dicMI, roundVal=self.roundVal)
            filename = filename.replace('_se.txt','') + '_summary.txt'
            if not self.hmi.entFile.write_file(self.rootTable, filename, striSummary, True):
                return False

        if self.showGraph or self.saveGraph:
            self.jsd.plot_JSD_MI(self, dicMI, self.filename_correction, xlabel, ylabel, onlyRef=False )

           
        return True

    def looping_JSD_one_sample(self):
        plt.close("all")
        plt.clf()

        self.jsd = JSD.Class_JSD(self)

        ''' -----------------------------------------'''  
        self.filename_correction = '_bias_corr' if self.withCorrection else ''
            
        self.time = crono.Timing()
        self.entFile = bioClass.Entropy_File()
        self.set_my_cluster_filenames()
        
        name = '%s %s %s'%(self.organism, self.gene, self.title)
        self.names = [[self.cluster_input_filename, name]]

        self.hmi = MI.Mutual_Information_Horizontal(self, want_nuc_aa_list=False, want_shannon = False)
        if not self.hmi.ok:
            self.error_msg = 'Problems with Init_MI_Heat_Map_Hist_Graph.'
            return False

        if self.organism == "Drosophila":
            self.dr = dro.Drosophila()
        
        # filename_dist_mat = self.cluster_input_filename
                   
        ''' always read and calc in looping_JSD_one_sample'''     
        # ret, dicMI, stri, sError  = self.read_and_calc_one_sample()
        ret, dicMI, _, _  = self.read_and_calc_one_sample()
       
        if not ret:
            return False
        
        xlabel = 'species x species'
        
        ylabel = 'JSD(%s)'%(self.unit)

        if self.showGraph or self.saveGraph:
            self.jsd.plot_JSD_MI(self, dicMI, self.filename_correction, xlabel, ylabel, onlyRef=True, stay=True)

        
        return True
    def set_my_cluster_filenames(self):
        self.set_gene_title()
        
        if self.withCorrection:
            self.filename_correction = '_bias_corr'
            self.str_correction = '-bias corr.'
        else:
            self.filename_correction = ''
            self.str_correction = ''
            

        if self.frame == 0:
            self.str_frame1 = ''
            self.str_frame2 = ''
        else:
            self.str_frame1 = 'frame=%i_'%(self.frame)
            self.str_frame2 = ' frame=%i, '%(self.frame)

        self.sufixTitle = '%s %s%s %s, gene=%s%s NOL=%i' %\
                (self.organism, self.minmax, self.str_correction, self.seqType, self.graph_gene_title, self.str_frame2, self.numOfLetters)
        self.sufix =     ('%s_%s_%s_%s_%sNOL%i_%iL_cutoff%i') %\
              (self.organism, self.minmax, self.seqType, self.gene_title, self.str_frame1, self.numOfLetters, self.cutoffLength, self.cutoffNumSeq)

        self.prefix = 'JSD_HMI_Sample_'
        
        if self.refSpecies:
            self.title_jsd = 'JSD(Horizontal Mutual Information) %s - %s\n%s'%\
                (self.unit, self.refSpecies, self.sufixTitle)
        else:
            if self.simulation: stri = "Simulation"
            elif self.oneSample: stri = "Sample"
            else: stri  = ""
            self.title_jsd = 'JSD(Horizontal Mutual Information) %s %s\n%s'%\
                (stri, self.unit, self.sufixTitle)
           

        self.cluster_input_filename = 'JSD_HMI_Sample_%s%s.txt' % (self.sufix, self.filename_correction)
        self.striMatrix = "Horizontal MI distance matrix"


    ''' looping_JSD_all_samples calls '''
    def read_and_calc_all_samples(self, showmessage=False):
        seqSpecies = []
        ''' stores a list of MI's for each species '''
        list_distSE = []
        list_distMI = []
        ns = []

        
        listSpecies = []
        listSpecJSD = []
        
        speciesList = sorted(self.all_species.keys())
        k = len(self.all_species.keys())

        str_sim = "sim_" if self.simulation else ""
        
        for iLoop in range(k) :
            if self.simulation:
                species = speciesList[iLoop]
                setSample = True
                #mat = self.all_species[species]
                # N = mat[0]
                #sampledSeqs = mat[1]
                #H0a = mat[2]
                sampleStr = "*"                
            else:
                species = speciesList[iLoop]
                mat = self.all_species[species]
                sampledSeqs = mat[1]
                setSample = True if sampledSeqs else False
                sampleStr = ""                

                            
            sp = self.which_sp(self.organism, species)

            numOfSeqs = self.sampleElems
            L = self.LH0

            if showmessage:
                stri = '%i/%i) %s #Seq=%i L=%i'%(iLoop+1, len(self.all_species.keys())+1, sp, numOfSeqs, L)
                self.showmsg_obs(stri, False)

            listSpecies.append(sp+sampleStr)
              

            if numOfSeqs == 0 or L == 0:
                self.error_msg = "numOfSeqs=%i, L=%i, for %s"%(numOfSeqs, L, species)
                continue
                          
            seqSpecies.append([species,sp,numOfSeqs])
            ns.append(numOfSeqs)

            self.hmi.ent.mySequence.set_seq_all_param(self, str_sim+species, seqs=None,setSample=setSample)
            self.hmi.ent.mySequence.lenSequence = L
            
            if showmessage:
                stri =  ' (%s) >> #seqs: %i, len=%i' % (sp, numOfSeqs, L)
                self.showmsg_obs(stri)

            sufix = self.hmi.ent.mySequence.sufix_hmi
            filenameMI = 'HMI_%s_mi.txt'%(sufix)
                            
            # print '>> %s sample len = %i ***** ' % (filenameMI, numSeqs)
            ''' hmi_Drosophila_americana_Gene_frame0_NOL1_400L_cutoff10_mij.txt '''
            if not self.entFile.setNames(self.rootTable, filename=filenameMI):
                self.showmsg_obs('Could not find %s'%(self.rootTable))
                return False, None, None, None

            if showmessage:
                self.showmsg_obs('Reading %s'%(filenameMI) )
                _ = self.time.start()

            ok, msg,  ValList, ValCorrList, _ = \
            self.entFile.read_HMI_files(root=self.rootTable, sufix=sufix, which='mi' ,showKeys=self.showmessage)

            if not ok:
                self.showmsg_obs(msg)
                return False, None, None, None

            ok, msg,  SElist, SEcorrList, _ = \
            self.entFile.read_HMI_files(root=self.rootTable, sufix=sufix, which='se', showKeys=self.showmessage)

            if not ok:
                self.showmsg_obs(msg)
                return False, None, None, None
               
                                     
            if showmessage:                       
                _ = self.time.finish()
                print self.time.milli(), 'ms'
            
            if not self.withCorrection:           
                arrayVal = np.array(ValList)
                arraySE = np.array(SElist)
            else:
                arrayVal = np.array(ValCorrList)
                arraySE = np.array(SEcorrList)


            sp = self.which_sp(self.organism, species)

            listSpecJSD.append(sp)
               
            ''' NORMALIZATION '''
            arrayVal /= np.sum(arrayVal)
            arraySE *= arrayVal
                                                
            list_distMI.append(arrayVal)
            list_distSE.append(arraySE)
            
            iLoop += 1

        dicMI, striJSD, sError = self.jsd.calc_JSD(listSpecies, list_distMI, list_distSE, parSE=self.numOfSEs,unit=self.unit)
        
        return True, dicMI, striJSD, sError
    
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
            elif species == self.refSpecies:
                stri = '%i/%i) %s #Seq=%i complement of %s'%(iLoop+1, len(self.all_species.keys())+1, species, N, str(sampledSeqs))
            else:
                stri = '%i/%i) %s #Seq=%i all sequences'%(iLoop+1, len(self.all_species.keys())+1, species, N)

            print stri
            
            sp = self.which_sp(self.organism, species)

            listSpecJSD.append(sp+sampleStr)
            
            numOfSeqs = self.sampleElems
            L = self.LH0

            listSpecies.append(sp+sampleStr)

            if numOfSeqs == 0 or L == 0:
                self.error_msg = "numOfSeqs=%i, L=%i, for %s"%(numOfSeqs, L, species)
                continue
                          
            seqSpecies.append([species,sp,numOfSeqs])
            ns.append(numOfSeqs)

            self.hmi.ent.mySequence.set_seq_all_param(self, species, seqs=None,setSample=setSample)
            self.hmi.ent.mySequence.lenSequence = L
            
            sufix = self.hmi.ent.mySequence.sufix_hmi
            filenameMI = 'HMI_%s_mi.txt'%(sufix) # sampleStr
                            
            # print '>> %s sample len = %i ***** ' % (filenameMI, numSeqs)
            ''' hmi_Drosophila_americana_Gene_frame0_NOL1_400L_cutoff10_mij.txt '''
            if not self.entFile.setNames(self.rootTable, filename=filenameMI):
                self.showmsg_obs('Could not find %s'%(self.rootTable))
                return False, None, None, None

            # _ = self.time.start()

            if (iLoop == 0 or species == self.refSpecies) and self.oneSample:
                ValList, ValCorrList, SElist, SEcorrList, _ = self.redo_hmi(species, isReference, N, sampledSeqs, hypothesis)
            else:
                ok, _,  ValList, ValCorrList, _ = \
                self.entFile.read_HMI_files(root=self.rootTable, sufix=sufix, which='mi' ,showKeys=self.showmessage)
    
                if not ok:
                    print "recalc: redo_hmi"
                    ValList, ValCorrList, SElist, SEcorrList, _ = self.redo_hmi(species, isReference, N, sampledSeqs, hypothesis)
                else:
                    ok, _,  SElist, SEcorrList, _ = \
                    self.entFile.read_HMI_files(root=self.rootTable, sufix=sufix, which='se', showKeys=self.showmessage)
        
                    if not ok:
                        print "recalc: redo_hmi"
                        ValList, ValCorrList, SElist, SEcorrList, _ = self.redo_hmi(species, isReference, N, sampledSeqs, hypothesis)
           
                    
            #_ = self.time.finish()
            #print self.time.milli(), 'ms'
        
            if not self.withCorrection:           
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


        dicMI, striJSD, sError = self.jsd.calc_JSD(listSpecies, list_distMI, list_distSE, parSE=self.numOfSEs,unit=self.unit)


        self.result = {}
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
             
        
        mean = np.mean(vals)*self.factor
        self.result["mean"] = round(mean,1)
        
        pvalue_inf, median = mquantiles(vals, prob=[0.05, 0.5])*self.factor
        self.result["pvalue_inf"] = round(pvalue_inf,1)
        self.result["median"] = round(median,1) 

        # print "\n\nSelected sample distribution -------------------"
        # print "mean=%.1f, median=%.1f, pvalue_inf=%.1f"%(self.result["mean"], self.result["median"], self.result["pvalue_inf"])
        # print "------------------------------------------------"
        
        if d == -1:
            stri = "Nothing found."
        else:
            d = round(d*self.factor,1)
            
            if self.refSpeciesH0a not in self.contig_table.keys():
                self.contig_table[self.refSpeciesH0a] = [0, 0]
                
            if d <= pvalue_inf:
                if self.refSpeciesH0a == "H0":
                    self.contig_table[self.refSpeciesH0a][0] += 1
                else:
                    self.contig_table[self.refSpeciesH0a][1] += 1
                stri = "%s is statistically significant with d=%.1f and p-value=%.1f"%(strRef, d, self.result["pvalue_inf"])
            else:
                if self.refSpeciesH0a == "H0":
                    self.contig_table[self.refSpeciesH0a][1] += 1
                else:
                    self.contig_table[self.refSpeciesH0a][0] += 1
                stri = "%s is not statistically significant with d=%.1f and p-value=%.1f"%(strRef, d, self.result["pvalue_inf"])

        print stri
        self.result["result"] = stri
        
        return True, dicMI, striJSD, sError


    def read_distance_matrix(self, filename):
        ret, _, colHeaders, dataMatrix = self.open_distance_matrix_file(self.rootTable + filename)
        if not ret:
            self.error_msg = 'Could not find %s'%(self.rootTable + filename)
            return False, None
        
        filename = filename.replace('.txt','') + '_se.txt'
        ret, _,_, seMatrix =  self.open_distance_matrix_file(self.rootTable + filename)
        if not ret:
            self.error_msg = 'Could not find %s'%(self.rootTable + filename)
            return False, None
        
                
        listSpecies = colHeaders
        numSpecies = len(colHeaders)

        dic = {}
        stri = 'species(%s)'%(self.unit)
        sError = 'species'
                
        for i in range(numSpecies):
            stri += '\t' + listSpecies[i]
            sError += '\t' + listSpecies[i]
                            
        for i in range(numSpecies):
            stri += '\n' + listSpecies[i]
            sError += '\n' + listSpecies[i]

            for j in range(numSpecies):
                if i == j:
                    dist, dError = 0.0, 0.0
                else:
                    spec1 = listSpecies[i] 
                    spec2 = listSpecies[j] 
                    
                    if i > j:
                        dist = dataMatrix[j][i]
                        dError = seMatrix[j][i]
                        
                        dic[spec2+' x '+spec1] = [dist,dError]
                    else:
                        dist = dataMatrix[i][j]
                        dError = seMatrix[i][j]
                    
                        dic[spec1+' x '+spec2] = [dist,dError]
        
                stri += '\t' + str(dist)
                sError += '\t' + str(dError)

        
        return True, dic, stri, sError
            
            
    def open_distance_matrix_file(self, filename):
            #open the file assuming the data above is in a file called 'dataFile'
            if not os.path.exists(filename):
                return False, None, None, None
                
            # print 'opening %s'%(filename)
        
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
            self.N = selSpec[0]
            self.sampledSeqs = selSpec[1]
            self.refSpeciesH0a = selSpec[2]
            
            self.spRef = self.which_sp(self.organism, spec)

            # self.showGraph = True if i <= 3 else False

            self.looping_JSD_one_sample()
            i += 1

        
                   
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
                
        return self.contig_table

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
    
    def save_dic(self, filename, dics):
        self.entFile.write_file_data_dic2(filename, dics)
        