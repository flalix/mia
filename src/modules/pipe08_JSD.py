#-*- coding: utf-8 -*-
'''
Created on 18/08/2014
Updated on 27/08/2014
Updated on 29/09/2014 -- Calc JSD by diff distribution or error propagation (not reviewed)

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica
'''
#import FileDialog
import classes.Mutual_Information_Vertical as imi
import classes.BioPythonClass as bioClass
import classes.Timing as timePack
import os
import classes.JSD as JSD
import numpy as np
from Tkinter import END
import matplotlib.pyplot as plt

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
        
        if desk.gene == '' and desk.title == '':
            self.error_msg = 'Define at least one Gene or Title'
            return

        self.jsd = JSD.Class_JSD(desk)

        if desk.mnat:
            desk.unit = 'mnat'
            desk.factor = 1000
            desk.roundVal = 4
        else:
            desk.unit = 'nat'        
            desk.factor = 1
            desk.roundVal = 6
            
        minmax_backup = desk.minmax
        correc_backup = desk.withCorrection
        
        if desk.each_all == 'each':
            matLoopList = [ [desk.minmax]]
        else:
            matLoopList = [ [kMinMax] for kMinMax in ['mincut','maxmer']]

    
        for matLoop in matLoopList:
            if not self.looping(desk, matLoop):
                desk.minmax = minmax_backup
                desk.withCorrection = correc_backup
                return

        desk.minmax = minmax_backup
        desk.withCorrection = correc_backup
        
        self.failure = False
        self.error_msg = 'Task ended. All right.'
        return

    def looping(self, desk, opt):
        print("-->>", opt)
        plt.close("all")
        plt.clf()

        desk.minmax =  opt[0]

        if desk.each_all == 'each':
            loopCorrection = [desk.withCorrection]
        else:
            loopCorrection = [True, False]            
        
        for correction in loopCorrection:
            ''' -----------------------------------------'''  
            desk.withCorrection = correction
                
            self.time = timePack.Timing()
            
            name = '%s %s %s'%(desk.organism, desk.gene, desk.title)
            print 'name', name
    
            self.entFile = bioClass.Entropy_File()
    
            desk.set_cluster_filenames() 
            
            desk.names = [desk.cluster_input_filename, name]
    
            self.imi = imi.Mutual_Information_Vertical(desk)
            if not self.imi.ok:
                self.error_msg = 'Problems with Init_MI_Heat_Map_Hist_Graph.'
                return False
            
            filename_dist_mat = desk.cluster_input_filename
                            
            if not os.path.exists(desk.rootTable +  filename_dist_mat) or desk.de_novo:
                ret, dicMI, stri, sError  = self.read_and_calc(desk)
            else:
                ret, dicMI, stri, sError = self.read_distance_matrix(desk, filename_dist_mat)
            
            if not ret:
                return False
            
            xlabel = 'species x species'

            ''' recalc: desk.title_jsd and desk.unit '''
            desk.set_cluster_filenames()
        
            ylabel = 'JSD(%s)'%(desk.unit)
                        
            desk.showmsg_obs('\n ============ %s ============'%(desk.striMatrix))
            desk.showmsg_obs(stri)
            desk.showmsg_obs('============ SE %s ============'%(desk.striMatrix))
            desk.showmsg_obs(sError)
            desk.showmsg_obs('============ %s ============\n'%(desk.striMatrix))

   
            if desk.saveData:
                if not self.imi.entFile.write_file(desk.rootTable, filename_dist_mat, stri, True):
                    return False
                
                filename = filename_dist_mat.replace('.txt','') + '_se.txt'
                if not self.imi.entFile.write_file(desk.rootTable, filename, sError):
                    return False
    
                striSummary = self.imi.print_MI_species_data(dicMI, roundVal=desk.roundVal)
                filename = filename.replace('_se.txt','') + '_summary.txt'
                if not self.imi.entFile.write_file(desk.rootTable, filename, striSummary, True):
                    return False
    
    
            listSpecies = self.get_species(desk)
            self.imi.ttest_MI_species_data(desk, listSpecies, dicMI)
            
            if desk.showGraph or desk.saveGraph:
                self.jsd.plot_JSD_MI(desk, dicMI, xlabel, ylabel)

           
        return True

    def get_species(self, desk):
        lista = desk.speciesListbox.get(0,END)
        listSpecies = []
       
        for line in lista:
            species = desk.find_species(line)
            
            mat = desk.dicParams[species]

            if mat[0] == 'x':
                listSpecies.append(species)
                
        return listSpecies
    

    def read_and_calc(self, desk):
        seqSpecies = []
        ''' stores a list of MI's for each species '''
        list_distSE = []
        list_distMI = []
        ns = []
                
        lista = desk.speciesListbox.get(0,END)

        iLoop = 0
        listSpecies = []
        listSpecJSD = []
        
        for line in lista:
            iLoop += 1
            species = desk.find_species(line)
            
            if desk.dr:
                sp = desk.dr.mnemonic(species)
            else:
                lista = species.split("_")
                sp = ''
                lim = 4
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

                        
            mat = desk.dicParams[species]

            if mat[0] != 'x':
                desk.showmsg_obs('%i/%i - %s not chosen - not checked)'%(iLoop, len(lista), species))
                continue

            if desk.minmax == 'mincut':
                numOfSeqs = mat[2]
                L = mat[5]
            else:
                numOfSeqs = mat[3]
                L = mat[6]

                
            stri = '%i/%i) %s #Seq=%i L=%i'%(iLoop, len(lista), sp, numOfSeqs, L)
            desk.showmsg_obs(stri, False)

            listSpecies.append(sp)
              

            if numOfSeqs == 0 or L == 0:
                self.error_msg = "numOfSeqs=%i, L=%i, for %s"%(numOfSeqs, L, species)
                continue
                          
            seqSpecies.append([species,sp,numOfSeqs])
            ns.append(numOfSeqs)
      
            
            self.imi.ent.mySequence.set_seq_all_param(desk, species, seqs=None, which=desk.which)
            self.imi.ent.mySequence.lenSequence = L
            stri =  ' (%s) >> #seqs: %i, len=%i' % (sp, numOfSeqs, L)
            desk.showmsg_obs(stri)

            if desk.which ==  'VMI':
                sufix = self.imi.ent.mySequence.sufix_vmi
                filenameMI = 'VMI_%s_mi.txt'%(sufix)
            elif desk.which ==  'VSH':
                sufix = self.imi.ent.mySequence.sufix_vmi
                filenameMI = 'VSH_%s_hShannon.txt'%(sufix)
            else:
                sufix = self.imi.ent.mySequence.sufix_hmi
                filenameMI = 'HMI_%s_mi.txt'%(sufix)
                                
            # print '>> %s sample len = %i ***** ' % (filenameMI, numSeqs)
            ''' hmi_Drosophila_americana_Gene_frame0_NOL1_400L_cutoff10_mij.txt '''
            if not self.entFile.setNames(desk.rootTable, filename=filenameMI):
                desk.showmsg_obs('Could not find %s'%(desk.rootTable))
                return False, None, None, None

            desk.showmsg_obs('Reading %s'%(filenameMI) )

            _ = self.time.start()
            
            if desk.which ==  'VMI':  # vertical mutual information
               
                ok, msg, ValList, ValCorrList = \
                    self.entFile.read_VMI_files(desk.rootTable, sufix, 'mi', showmessage=False)

                if not ok:
                    desk.showmsg_obs(msg)
                    return False, None, None, None
                             
                ok, msg, SElist, SEcorrList = \
                    self.entFile.read_VMI_files(desk.rootTable, sufix, 'se', showmessage=False)

                if not ok:
                    desk.showmsg_obs(msg)
                    return False, None, None, None
            elif desk.which ==  'VSH':  # vertical shannon entropy
               
                ok, msg, ValList, ValCorrList = \
                    self.entFile.read_VMI_files(desk.rootTable, sufix, 'hShannon', showmessage=False)

                if not ok:
                    desk.showmsg_obs(msg)
                    return False, None, None, None
                             
                ok, msg, SElist, SEcorrList = \
                    self.entFile.read_VMI_files(desk.rootTable, sufix, 'seh', showmessage=False)

                if not ok:
                    desk.showmsg_obs(msg)
                    return False, None, None, None
            else:  # horizontal mutual information

                ok, msg,  ValList, ValCorrList, _ = \
                self.entFile.read_HMI_files(root=desk.rootTable, sufix=sufix, which='mi' ,showKeys=desk.showmessage)

                if not ok:
                    desk.showmsg_obs(msg)
                    return False, None, None, None

                ok, msg,  SElist, SEcorrList, _ = \
                self.entFile.read_HMI_files(root=desk.rootTable, sufix=sufix, which='se', showKeys=desk.showmessage)

                if not ok:
                    desk.showmsg_obs(msg)
                    return False, None, None, None
                                                            
            _ = self.time.finish()
            print self.time.milli(), 'ms'
            
            if not desk.withCorrection:           
                arrayVal = np.array(ValList)
                arraySE = np.array(SElist)
            else:
                arrayVal = np.array(ValCorrList)
                arraySE = np.array(SEcorrList)


            if desk.dr:
                sp = desk.dr.mnemonic(species.replace('Drosophila ',''))
            else:
                sp = species

            listSpecJSD.append(sp)
               
            ''' NORMALIZATION '''
            arrayVal /= np.sum(arrayVal)
            arraySE *= arrayVal
                                                
            list_distMI.append(arrayVal)
            list_distSE.append(arraySE)

        dicMI, stri, sError = self.jsd.calc_JSD(listSpecies, list_distMI, list_distSE)

        return True, dicMI, stri, sError
    
    
    def read_distance_matrix(self, desk, filename):
        ret, _, colHeaders, dataMatrix = self.open_distance_matrix_file(desk.rootTable + filename)
        if not ret:
            self.error_msg = 'Could not find %s'%(desk.rootTable + filename)
            return False, None
        
        filename = filename.replace('.txt','') + '_se.txt'
        ret, _,_, seMatrix =  self.open_distance_matrix_file(desk.rootTable + filename)
        if not ret:
            self.error_msg = 'Could not find %s'%(desk.rootTable + filename)
            return False, None
        
                
        listSpecies = colHeaders
        numSpecies = len(colHeaders)

        dic = {}
        stri = 'species(%s)'%(desk.unit)
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
        