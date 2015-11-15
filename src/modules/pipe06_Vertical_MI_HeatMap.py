'''
Created on 05/08/2013
Update  on 30/09/2013
Update  on 14/02/2014
Update  on 24/04/2014
Update  on 24/04/2014
Updated on 07/08/2014 - Pipe

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica
'''
#import FileDialog
import classes.Mutual_Information_Vertical as VMI
import classes.BarGraphic as graphPac
import numpy as np
import classes.Drosophila as dro
from Tkinter import END
import matplotlib.pyplot as plt

class Pipe():
    def __init__(self, desk, which_module):

        self.desk = desk
        self.failure = True
        self.error_msg = ''
        
        desk.which = 'VMI'
        
        
        try:
            desk.get_params()
        except:  
            self.error_msg = 'Could not get parameters.'
            return

                             
        ''' -----------------------------------------'''
        if desk.organism == '':
            self.error_msg = "Define the organism or write 'any'"
            return
        
        if desk.gene_title == '':
            self.error_msg = 'Define at least one Gene or Title'
            return
        
        desk.isProtein = not(desk.dna_prot == 'DNA') 
                    
        ''' -------------------'''
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
            if not self.looping(desk, matLoop, lista, which_module):
                return

        self.failure = False
        self.error_msg = 'Task ended. All right.'
        return

    def looping(self, desk, opt, lista, which_module):
        print("\n---------------------------------")
        plt.close("all")
        plt.clf()
                    
        desk.minmax =  opt[0]
        self.list_species_name_vmi = []

        iLoop = 0
        
        for line in lista:
            ''' Vertical_MI initialize without nuc_aa_List '''
            vmi = VMI.Mutual_Information_Vertical(desk,  want_nuc_aa_list=False)   
        
            if not vmi.ok:
                self.error_msg = 'Problems with Vertical_MI.'
                return False                
            
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
            
            filename = ('%s_%s_%s_%s_%s_%iL_cutoff%i_consensus.fasta') %\
                       (desk.organism, desk.minmax, species, desk.seqType, desk.gene_title, desk.cutoffLength, desk.cutoffNumSeq)

            name = ('%s_%s') %  (desk.organism, species)

            if numOfSeqs == 0 or L == 0:
                self.error_msg = "numOfSeqs=%i, L=%i, for %s"%(numOfSeqs, L, name)
                continue

            if not vmi.ent.mySequence.read_fasta_3methods(desk, filename, name, method=desk.label_random, showmessage=desk.showmessage):
                self.error_msg = "Could not find %s"%(filename)
                return False
        
            if numOfSeqs != vmi.ent.mySequence.numSequences:
                print 'Strange: #seqs disagree: saved %i x read %i'%(numOfSeqs, vmi.ent.mySequence.numSequences)
                return False
            
            if L != vmi.ent.mySequence.lenSequence:
                print 'Strange: length of seqs disagree: saved %i x read %i'%(L, vmi.ent.mySequence.lenSequence)
                return False


            seqs = []
            for rec in vmi.ent.mySequence.seqs:
                seqs.append(str(rec)) 
             
            '''  all params in mySequence '''
            vmi.ent.mySequence.set_seq_all_param(desk, species, seqs, which="VMI")

            sufix = vmi.ent.mySequence.sufix_vmi
            mat = [filename, name, species, vmi, sufix]
            self.list_species_name_vmi.append(mat)
            
            if desk.each_all == 'each':
                ret, msg = vmi.prepareVerticalMutualInfo(desk, sufix=sufix, which_module=which_module, showmessage=desk.showmessage)
            else:
                ret, msg = vmi.prepareVerticalMutualInfo(desk, sufix=sufix, which_module='all', showmessage=desk.showmessage)
                
            if not ret:
                desk.showmsg_obs(msg)
                return False

        ''' end loop of selected species '''


        self.save_tables_show_graph(desk, which_module) 

        del self.list_species_name_vmi
        
        return True
    
    def save_tables_show_graph(self, desk, which_module):

        if desk.each_all == 'each':
            list_modules = [which_module]
            bias_corr_list = [desk.withCorrection]
        else:
            list_modules = ['Entropy', 'VMI']
            bias_corr_list = [False, True]
                    
        print("-----------------------------------------------------------")
        for which_module in list_modules:
            for bias_corr in bias_corr_list:
                iLoop = 0
                listMI_Anova = []
                speciesParams = []
                
                for mat in self.list_species_name_vmi:
                    iLoop += 1
                    filename, name, species, vmi, sufix = mat
                    
                    mat = desk.dicParams[species]
        
                    if desk.minmax == 'mincut':
                        numOfSeqs = mat[2]
                        L = mat[5]
                    else:
                        numOfSeqs = mat[3]
                        L = mat[6]
       
                    print("%i) %s - %s %s %s"%(iLoop, name, desk.minmax, "bias corr." if bias_corr else 'no corr.', which_module ))
        
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
    
                    std_rand = "" if (desk.label_random == "") else " " + desk.label_random
                    
                    if which_module == 'Entropy':
                        title = 'Entropy Distribution %s %s, %s %s\n%s%s%s; %i seqs; len=%i; #letter=%i'%\
                            (desk.organism, species, desk.seqType, desk.gene_title, desk.minmax, str_correction, std_rand, numOfSeqs, L,  desk.numOfLetters)
                    else:
                        title = 'VMI Heat Map %s %s, %s %s %s%s%s\n%i seqs; len=%i; #letter=%i'%\
                            (desk.organism, species, desk.seqType, desk.gene_title, desk.minmax, str_correction, std_rand, numOfSeqs, L,  desk.numOfLetters)
            
            
                    if which_module == 'Entropy':
                        if not desk.withCorrection:           
                            arrayVal = np.array(vmi.HShannonList)
                            arraySE = np.array(vmi.SeHShannonList)
                        else:
                            arrayVal = np.array(vmi.HShannonCorrList)
                            arraySE = np.array(vmi.SeHShannonCorrList)
                    else:
                        if not desk.withCorrection:           
                            arrayVal = np.array(vmi.MIlist)
                            arraySE = np.array(vmi.SeMIList)
                        else:
                            arrayVal = np.array(vmi.MIcorrList)
                            arraySE = np.array(vmi.SeMICorrList)
        
                    ''' normalization dividing by numOfLetters '''
                    if desk.norm:
                        arrayVal /= desk.numOfLetters
                        arraySE /= desk.numOfLetters
        
                    ''' mili nats '''
                    if desk.mnat:
                        arrayVal *= 1000
                        arraySE *= 1000
                        
                        roundVal = 2
                    else:
                        roundVal = 4
        
                    # print '--- params -----------------'
                    maxMI = 0
                    maxiPos = None
                    
                    is_zero = True
                    
                    for pos in range(len(arrayVal)):
                        if arrayVal[pos] > maxMI:
                            maxMI = arrayVal[pos]
                            SE = arraySE[pos]
                            i,j = vmi.ijList[pos]
                            maxiPos = [i,j,maxMI,SE]
                            is_zero = False
        
                    if is_zero:
                        stri = '### Species %s has MI = ZERO. Too conserved data sequences. Impossible to include in analysis.'%(species)
                        print(stri)
                        continue
                        
                    xSeq = [x for x in range(len(arrayVal))]
        
                    self.meanY = np.round(np.mean(arrayVal), roundVal)
                    self.medianY = np.round(np.median(arrayVal), roundVal)
                    self.stdY  = np.round(np.std(arrayVal), roundVal)
                    self.maxY = np.round(np.max(arrayVal), roundVal)
                    self.minY = np.round(np.min(arrayVal), roundVal)
        
                    speciesParams.append([species,sp,numOfSeqs, self.meanY, self.medianY, self.stdY, self.maxY, self.minY])
                    
                    if desk.saveData:
                        listMI_Anova.append(arrayVal)
        
                    if desk.showGraph or desk.saveGraph:
                        if which_module == 'Entropy':
                            gHist = graphPac.Histogram_FreqDistribution(desk, title)
                            
                            gHist.meanY = self.meanY
                            gHist.medianY = self.medianY
                            gHist.stdY  = self.stdY
                            gHist.maxY = self.maxY
                            gHist.minY = self.minY
                                            
                            gHist.plot_H_MI('V', xSeq=xSeq, ySeq=arrayVal, ySE=arraySE, showError=True)
                            gHist.densityBar('V', seq=arrayVal)
                            
                            pictureName = 'VHShannon_%s%s%s'%(sufix, filename_correction, desk.stri_random)
                            gHist.myPlot.print_graph(desk, gHist.fig, pictureName, frame=desk.tk_root)
                        else:
                            ''' 3D dont has ceil '''
                            if desk.is3D:
                                limSup = self.maxY
                            else:
                                ceil = desk.heatmap_ceil_value
                                
                                ''' the same roof for all heatmaps: parametrize in future '''
                                if desk.heatmap_ceil:
                                    if self.maxY <= ceil:
                                        limSup = ceil
                                    else:
                                        limSup = self.maxY
                                else:
                                    limSup = self.maxY
                            
                            ''''  updated 28/09/2015 '''
                            gMI = vmi.plotHeatMap(desk, desk.is3D, arrayVal, vmi.ijList, L, maxiPos, title, species=species, limSup=limSup, roundVal=roundVal, str_correction=str_correction)
                 
                            if desk.is3D:
                                pictureName = 'HeatMap_3D_VMI_%s%s%s'%(sufix, filename_correction,desk.stri_random)
                            else:
                                gMI.densityHeatmapBar(arrayVal, limSup)
                                pictureName = 'HeatMap_2D_VMI_%s%s%s'%(sufix, filename_correction,desk.stri_random)

                            gMI.myPlot.print_graph(desk, gMI.fig, pictureName, frame=desk.tk_root)
        
                            '''
                            plt.cla()
                            plt.clf()
                            plt.close()
                            del gMI
                            gc.collect()
                            '''

                # Flavio 02/06/2015
                if which_module == 'Entropy':
                    which_symb = 'VHShannon'
                else:
                    which_symb = 'VMI'
                    
                    
                if desk.saveData:
                    sufix = ('%s_%s_%s_%s_NOL%i_%iL_cutoff%i') %\
                            (desk.organism, desk.minmax, desk.seqType, desk.gene_title, desk.numOfLetters, desk.cutoffLength, desk.cutoffNumSeq)
            
                    filename = '%s_summary_%s%s%s.txt' % (which_symb, sufix, filename_correction, desk.stri_random)
                    stri = vmi.ent.calc_anova(listMI_Anova, sufix)
                    # desk, speciesParams, roundVal=4, filename=None, stri = '', saveData=False):
                    vmi.ent.print_data_summary(desk, speciesParams, roundVal=roundVal, filename=filename, stri=stri)
                            