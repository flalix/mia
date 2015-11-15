#-*- coding: utf-8 -*-
'''
Created on 16/04/2013
Updated on 25/04/2013
Updated on 14/10/2013
Updated on 30/01/2014
Updated on 24/04/2014
Updated on 13/05/2014
Updated on 02/06/2014
Updated on 16/06/2014 - with Consensus
Updated on 28/07/2014 - Pipe
Updated on 09/10/2014 - new MI_Result with correction
Updated on 27/09/2015 - shuffling or randomizing

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica
'''
# import FileDialog
import classes.Mutual_Information_Horizontal as MI
import classes.BarGraphic as graphPac
import classes.Drosophila as dro
from Tkinter import END
import numpy as np
import matplotlib.pyplot as plt
import gc
from numpy.random.mtrand import shuffle

class Pipe():
    def __init__(self, desk):

        self.desk = desk
        self.failure = True
        self.error_msg = ''
        
        try:
            desk.get_params()
        except:  
            self.error_msg = 'Could not get parameters.'
            return
        
        
        if desk.organism == '':
            self.error_msg = "Define the organism or write 'any'"
            return
        
        if desk.gene_title == '':
            self.error_msg = 'Define at least one Gene or Title'
            return

        ''' -----------------------------------------'''            
        
        ''' can be 0 any k, 1 for desk.frame 1, a
            nd > 1 displaying all frames and 0 '''
        
        self.arrColor = ['blue','Black','red','green','magenta','cyan','Yellow']
        self.arrLinestyleCode = ['k--','--','--','k--']


        if desk.organism == "Drosophila":
            self.dr = dro.Drosophila()

                 

        lista = desk.speciesListbox.get(0,END)
        
        ''' must choose at least one species '''
        at_least_one = False
        
        for line in lista:
            species = desk.find_species(line)
            mat = desk.dicParams[species]

            if mat[0] == 'x':
                at_least_one = True
                break
                                 
        if not at_least_one:
            self.failure = False
            desk.showmsg_obs('Please select at least one sequence.')
            return


        ''' ------------------------------------------------------- '''      
                 
        if desk.each_all == 'each':
            matLoopList = [ [desk.minmax]]
        else:
            matLoopList = [ [kMinMax] for kMinMax in ['mincut','maxmer'] ]

    
        for matLoop in matLoopList:
            if not self.looping(desk, matLoop, lista):
                return
        
        self.failure = False
        self.error_msg = 'Task ended. All right.'
        return

    def looping(self, desk, opt, lista):
        print("\n---------------------------------")
        plt.close("all")
        plt.clf()
                    
        desk.minmax =  opt[0]
        self.list_species_name_hmi = []
        
        
        iLoop = 0
        
        for line in lista:
            hmi = MI.Mutual_Information_Horizontal(desk, want_nuc_aa_list=False)
            
            if not hmi.ok:
                return
            
            iLoop += 1
            species = desk.find_species(line)
            mat = desk.dicParams[species]

            if mat[0] != 'x':
                desk.showmsg_obs('%i/%i - %s not chosen - not checked.'%(iLoop, len(lista), species))
                continue

            if desk.minmax == 'mincut':
                numOfSeqs = mat[2]
                L = mat[5]
            else:
                numOfSeqs = mat[3]
                L = mat[6]
                
            stri = '%i/%i) %s #Seq=%i L=%i'%(iLoop, len(lista), species, numOfSeqs, L)
            desk.showmsg_obs(stri, True)

            ''' reading each consensus file '''
            filename = ('%s_%s_%s_%s_%s_%iL_cutoff%i_consensus.fasta') %\
                 (desk.organism, desk.minmax, species, desk.seqType, desk.gene_title, desk.cutoffLength, desk.cutoffNumSeq)

            name = '%s_%s%s'%(desk.organism, species, desk.gene_title)

            if numOfSeqs == 0 or L == 0:
                self.error_msg = "numOfSeqs=%i, L=%i, for %s"%(numOfSeqs, L, name)
                continue
            
            ''' shuffling or randomizing '''
            if not hmi.ent.mySequence.read_fasta_3methods(desk, filename, name, method=desk.label_random, showmessage=desk.showmessage):
                self.error_msg = "Could not find %s"%(filename)
                return False


            if numOfSeqs != hmi.ent.mySequence.numSequences:
                print 'Strange: #seqs disagree: saved %i x read %i'%(numOfSeqs, hmi.ent.mySequence.numSequences)
                return False
            
            if L != hmi.ent.mySequence.lenSequence:
                print 'Strange: length of seqs disagree: saved %i x read %i'%(L, hmi.ent.mySequence.lenSequence)
                return False
        
            seqs = []
            for rec in hmi.ent.mySequence.seqs:
                seqs.append(str(rec))
                                    
            '''  all params in mySequence '''
            desk.which = 'HMI'
            hmi.ent.mySequence.set_seq_all_param(desk, species, seqs, which="HMI")

            # eFile = bpc.Entropy_File()
            
            sufix = hmi.ent.mySequence.sufix_hmi
            mat = [filename, name, species, hmi, sufix]
            self.list_species_name_hmi.append(mat)
            
            ret, msg = hmi.prepareHorizontalMutualInfo(desk, sufix=sufix, consensus=True, showmessage=desk.showmessage)
                
            if not ret:
                desk.showmsg_obs(msg)
                return False
            
            
        ''' end loop of selected species '''

        self.save_tables_show_graph(desk) 

        del self.list_species_name_hmi
        
        try:
            gc.collect() 
        except:
            pass    

                  
        return True              
            
            
    def save_tables_show_graph(self, desk):

        if desk.each_all == 'each':
            bias_corr_list = [desk.withCorrection]
        else:
            bias_corr_list = [False, True]

        print("-----------------------------------------------------------")
        for bias_corr in bias_corr_list:
            iLoop = 0
            speciesParams = []
            listMI_Anova = []       
                 
            for mat in self.list_species_name_hmi:
                iLoop += 1
                filename, name, species, hmi, sufix = mat
                
                mat = desk.dicParams[species]
    
                if desk.minmax == 'mincut':
                    numOfSeqs = mat[2]
                    L = mat[5]
                else:
                    numOfSeqs = mat[3]
                    L = mat[6]
   
                print("%i) %s - %s, %s"%(iLoop, name, desk.minmax, "bias corr." if bias_corr else 'no corr.' ))
    
                if desk.organism == 'Drosophila':
                    sp = self.dr.mnemonic(species.replace('Drosophila ',''))
                else:
                    sp = species
                        
                desk.withCorrection = bias_corr
                
                if desk.withCorrection:
                    str_correction = ' (bias corr.)'
                    filename_correction = '_bias_corr'
                else:
                    str_correction = ''
                    filename_correction = ''        

                if desk.frame == 0:
                    sFrame = ''
                else:
                    sFrame = ', frame=' + str(desk.frame)

    
                std_rand = "" if (desk.label_random == "") else " " + desk.label_random
                
                    
                title = 'HMI: %s %s, %s %s\n%s%s%s; %i seqs; len=%i%s; #letter=%i'%\
                        (desk.organism, species, desk.seqType, desk.gene_title, desk.minmax, str_correction, std_rand, numOfSeqs, L, sFrame, desk.numOfLetters)
    
                                
                ''' HMI mean and sdv for each k from n sequences '''
                if not desk.withCorrection:           
                    arrayMIthis = np.array(hmi.arrayMI)
                    arraySE     = np.array(hmi.arraySE)
                else:
                    arrayMIthis = np.array(hmi.arrayMIcorr)
                    arraySE     = np.array(hmi.arraySEcorr)
    
               
                if desk.mnat:
                    arrayMIthis *= 1000
                    arraySE *= 1000
    
                    roundVal = 1
                    unit = 'mnat'
                else:
                    roundVal = 4
                    unit = 'nat'
                    
                if desk.norm:
                    arrayMIthis /= desk.numOfLetters
                    arraySE /= desk.numOfLetters
                
                # desk.showmsg_obs('  >>> species %s with %i sequences'%(species, numOfSeqs))
                '''  first time build abcissa
                         if choose == 1 and desk.frame == 1
                         or desk.frame == 0
                '''
                xSeq = []
                
                if desk.frame < 2:
                    for k in range( len(arrayMIthis)):
                        xSeq.append(k+3)
                else:
                    ''' frame1: k = {3,6,9, ... 3n} '''
                    for k in range( len(arrayMIthis)):
                        xSeq.append((k+1)*3)
                                        
           
                arrayVal = []
                ySeqSEFinal = []


                ''' if only one curve '''
                if desk.frame < 2:
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
                
                if desk.showGraph or desk.saveGraph:
                    if desk.frame < 2:
                        gHist = graphPac.Histogram_FreqDistribution(desk, title)
                        
                        gHist.meanY = self.meanY
                        gHist.medianY = self.medianY
                        gHist.stdY  = self.stdY
                        gHist.maxY = self.maxY
                        gHist.minY = self.minY
                                        
                        gHist.plot_H_MI('V', xSeq=xSeq, ySeq=arrayVal, ySE=arraySE, showError=True)                
                        
                        gHist.plot_H_MI('H', xSeq=xSeq, ySeq=arrayVal, ySE=ySeqSEFinal, showError=True)
                        
                        gHist.densityBar("H", seq=arrayVal)
                    else:
                        gHist.sameBar(xSeq=xSeq, arrayMIthis=arrayVal,linestyleCode=self.arrLinestyleCode[desk.frame], color=self.arrColor[desk.frame])
                            
                    
                    if (desk.frame == 0) or (desk.frame==3):    
                        pictureName = 'HMI_%s%s'%(sufix, filename_correction)
                        gHist.myPlot.print_graph(desk, gHist.fig, pictureName, frame=desk.tk_root)
    
                    '''
                    del gHist
                    plt.cla()
                    plt.clf()
                    gc.collect()
                    '''
                    
            if desk.saveData:
                sufix = ('%s_%s_%s_%s_frame%i_NOL%i_%iL_cutoff%i') %\
                        (desk.organism, desk.minmax, desk.seqType, desk.gene_title, desk.frame, desk.numOfLetters, desk.cutoffLength, desk.cutoffNumSeq)

                std_rand = "" if (desk.label_random == "") else desk.stri_random

                filename = 'HMI_summary_%s%s%s.txt' % (sufix, filename_correction, std_rand)        
                stri = hmi.ent.calc_anova(listMI_Anova, sufix)
        
                hmi.ent.print_data_summary(desk, speciesParams, roundVal=roundVal, filename=filename, stri=stri)

        self.failure = False
        self.error_msg = 'Task ended. All right.'
        return True

