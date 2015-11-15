'''
Revised on 20/02/2014

@author: Flavio Lichtenstein
@local: UNIFESP - Bioinformatica
'''
import FileDialog
import classes.BioPythonClass as bioClass
import classes.BarGraphic as bg
import random
import numpy as np

class Pipe():
    def __init__(self, desk):
        self.desk = desk
        
        ''' ------ params -----'''
        self.failure = True
        self.error_msg = ''
        self.showing_graph = False
        
        try:
            desk.saveData = desk.saveData_var.get()
            desk.numOfLetters = desk.numOfLetters_var.get()
            desk.numOfExperiments = desk.numOfExperiments_var.get()

            desk.showGraph = desk.showGraph_var.get() 
            desk.saveGraph = desk.saveGraph_var.get()
            
            self.showing_graph = desk.showGraph
                
            desk.dna_prot = desk.dna_prot_var.get()
            desk.dpi      = desk.dpi_var.get()
            
            desk.rootFasta = desk.rootFasta_var.get()
            desk.rootImage = desk.rootImage_var.get()
            desk.rootTable = desk.rootTable_var.get()
            desk.rootEntropy = desk.rootEntropy_var.get()
          
            desk.showmessage = desk.showmessage_var.get()
          
        except:
            self.error_msg = 'Could not get parameters.' 
            return
        

        if desk.dna_prot == "DNA":
            title = "DNA"
            nuc_aa = "AGTC"
        else:
            title = "Amino Acids"
            nuc_aa = "ACDEFGHIKLMNPQRSTVWY"
        title += " Random Shannon Entropy Simulation with " + str(desk.numOfExperiments) + " of experiments, word of " + str(desk.numOfLetters) + " letters"
        
        start_at = desk.sartAt_var.get()
        numIndiv = desk.lengthSim_var.get()
        numRealIndiv = numIndiv - start_at + 1

        if start_at < desk.numOfLetters*2:
            print 'Increase start_at, to low'
            
        if numRealIndiv < desk.numOfLetters*20:
            print 'Increase numIndiv, to low'
            
        if numRealIndiv < 20:
            print 'numRealIndiv < 20, to low'
            exit()
        
        left = 0.02
        top = 0.94
        
        
        desk.showmsg_obs('--------------------------')
        desk.showmsg_obs('from %i to %i = %i'%(start_at, numIndiv, numRealIndiv) )
        desk.showmsg_obs('--------------------------')
        
        seqX = []                        
        for ii in range(start_at,numIndiv+1,desk.numOfLetters):
            seqX.append(ii+start_at)
            
        i = 0
        
        numCols=3
        numLines=1
        
        legColumns="Random Study"

        ent = bioClass.Entropy(desk=desk)
        
        meanShannonArray = np.zeros((numRealIndiv , 1))
        sdvShannonArray = np.zeros((numRealIndiv , 1))
        
        
        for numOfExp in range(desk.numOfExperiments):
            desk.showmsg_obs('%i/%i' %(numOfExp+1, desk.numOfExperiments))
            ArrAvgShannon = []
            ArrSEShannon = []
            
            count = start_at
            for i in range(start_at, numIndiv+1, desk.numOfLetters):
                ent.mySequence.clearSequences()
                ent.numIndiv = count
        
                seqs = []
        
                for _ in range(i):
                    stri = ''
                    for _ in range(desk.numOfLetters):
                        stri += random.choice(nuc_aa)
                    seqs.append(stri)
                        
                ''' -------------------------------------------------------
                --- Calculating entropy
                ------------------------------------------------------- '''
        
                '''
                    generates:
                    self.dicPiList
                    self.entropyShannon
                    self.varEntropyShannon
                '''        
                ent.calcHSannon_bias_correction(seqs, desk.numOfLetters)                

                ArrAvgShannon.append(np.mean(ent.HShannonList))
                ArrSEShannon.append(np.mean(ent.SeHShannonList)) 
                
                count += 1
                    
            for i in range(len(ArrAvgShannon)):
                meanShannonArray[i] += (ArrAvgShannon[i] / desk.numOfExperiments)
                sdvShannonArray[i] += (ArrSEShannon[i] / desk.numOfExperiments)
        
                
        print 'lenX', len(seqX)
        
        print '\n--- meanShannonArray ----------'
        for i in range(len(ArrAvgShannon)):
            print i*desk.numOfLetters, ArrAvgShannon[i], '(', ArrSEShannon[i], ') \t%Mer=', round(100*ArrAvgShannon[i]/(2*desk.numOfLetters),3)
             
        seqMean = []
        seqSdv = []
        seqVC = []
        
        
        for ii in range( len(ArrAvgShannon)):
            if meanShannonArray[ii] > 0:
                VC = round(  (sdvShannonArray[ii] / meanShannonArray[ii]) * 100.0,3)
            else:
                VC = float('inf')
                
                
            seqMean.append([meanShannonArray[ii],0,0,]) 
            seqSdv.append([sdvShannonArray[ii],0,0,]) 
            seqVC.append([VC,0,0,]) 
        
        
        printLegend = False
        label=""
        
        xlabel="n (num of seqs)"
        
        subTitle= "<H> Shannon"
        ylabel="<h>"
        
        graph = bg.MultiLineWithTitle(numLines, numCols, left, top, legColumns, legendTitle='', title=title,dpi=desk.dpi)
        
        graph.multilineGeneralWithTitle(0, 1, seqX, seqMean, title=subTitle, label=label, printLegend=printLegend, xlabel=xlabel, ylabel=ylabel)
        
        subTitle= "SE(<h>)"
        ylabel="SE"
        graph.multilineGeneralWithTitle(0, 2, seqX, seqSdv, title=subTitle, label=label, printLegend=printLegend, xlabel=xlabel, ylabel=ylabel)
        
        subTitle= "Variational Coefficient"
        ylabel="VC%"
        graph.multilineGeneralWithTitle(0, 3, seqX, seqVC, title=subTitle, label=label, printLegend=printLegend, xlabel=xlabel, ylabel=ylabel)
        
        
        entFile = bioClass.Entropy_File()

        seqVal  = []
        seqSdv  = []
        
        
        print "i \t <h>shannon \t  dp \t VC%"
        for i in range(numRealIndiv):
            # print meanShannonArray[i][0]
            mean = meanShannonArray[i][0]
            sdv  = sdvShannonArray[i][0]
            
            # value = [mean, mean+sdv, mean-sdv]
                
            seqVal.append(mean)
            seqSdv.append(sdv)
                    
                     
            if mean > 0:
                print str(i+start_at) + '\t' + str(mean) + '\t' + str(sdv) + '\t' + str( round( (sdv/mean) * 100, 5))
            else:
                print str(i+start_at) + '\t' + str(mean) + '\t' + str(sdv) + '\t None'
                 
              
        if desk.saveData:                               
            filename = 'shannon_random_DNA_Letter%i_Exp%i_dic.txt'%(desk.numOfLetters, desk.numOfExperiments)
                                           
            if entFile.setNames(desk.rootEntropy, filename):
                entFile.save_shannon_file(numIndiv, desk.numOfExperiments, seqX, seqVal, seqSdv)
        
        pictureName = 'shannon_random_DNA_Letter%i_Exp%i'%(desk.numOfLetters, desk.numOfExperiments)
        graph.myPlot.print_graph(self.desk, graph.fig, pictureName, frame=self.desk.tk_root, stay=False)
           
        
        self.failure = False
        return
        
