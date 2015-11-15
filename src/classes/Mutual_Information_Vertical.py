#-*- coding: utf-8 -*-
'''
Created on 16/04/2013
Updated on 18/04/2013
Updated on 07/11/2013
Updated on 28/02/2014
Updated on 08/05/2014


@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica

'''
import BioPythonClass as bioClass
import BarGraphic as graphPac
import Timing as crono
import os, numpy as np
import matplotlib.pyplot as plt
import math


class Mutual_Information_Vertical:
    def __init__(self, desk, want_nuc_aa_list=False, want_shannon = False):
        self.ok = False

        if want_nuc_aa_list:
            self.build_dna_aa(desk)
            
        self.entFile = bioClass.Entropy_File()

        self.legColumns = 1
        self.legendTitle = '%s Species'%(desk.organism)
        self.title = 'Vertical MI'
        
        ''' numLines=1, numCols=1, 
                       legColumns=self.legColumns, 
                       legendTitle=self.legendTitle, title=self.title '''
        self.ent = bioClass.Entropy(desk)

        
        if want_shannon:
            filenameShannon = 'shannon_random_DNA_Letter%i_Exp100_dic.txt'
            self.filenameShannon1 = filenameShannon%(desk.numOfLetters)
            self.filenameShannon2 = filenameShannon%(desk.numOfLetters*2)
            
            ok, msg = self.ent.read_shannon_random_file(desk.rootEntropy, self.filenameShannon1, 1, showKeys = False)
            if not ok:
                desk.showmsg_obs(msg + ", didn't find Shannon dictionary for 1 letter.")
       
            ok, msg = self.ent.read_shannon_random_file(desk.rootEntropy, self.filenameShannon2, 2, showKeys = False)
            if not ok:
                desk.showmsg_obs(msg + ", didn't find Shannon dictionary for 2 letters.")
                return
    
        self.ok = True

    def noZeros(self, soma, param):
        maxi = 0
        zeros = 0
        iMaxi = 0
        
        for i in range(len(soma)):
            if soma[i] == 0:
                zeros += 1
            else:
                if soma[i] > maxi:
                    maxi = soma[i]
                    iMaxi  = i
        
        if zeros == 0:
            return soma
        
        ''' each zero incremente delta
            from maximum decrement delta '''
        for i in range(len(soma)):
            if soma[i] == 0:
                soma[i] = param
                soma[iMaxi] -= param
        
        return soma
        

    def prepareVerticalMutualInfo(self, desk, sufix, which_module, showmessage=False):
        self.desk = desk
        de_novo = desk.de_novo

        # desk.showmsg_obs('Prepare Mutual Information: num indiv: %i length: %i' % (self.ent.mySequence.getNumOfSeqs(),self.ent.mySequence.getLengthSequence()))
        
        ''' setNames define completeFilename +++ SPECIE !'''
        if not self.entFile.setPrefixSufix(root=desk.rootTable, prefix='VMI', sufix=sufix):
            return

        filename = self.entFile.MI_prefix_sufix_txt(self.entFile.prefix_sufix, 'mi')
        filename = desk.rootTable + filename

        time = crono.Timing()
        time.start()

        if not os.path.exists(filename) or de_novo:
            stri = "new: lines %i cols %i, %s"%(len(self.ent.mySequence.seqs), len(self.ent.mySequence.seqs[0]), filename)
            desk.showmsg_obs(stri)
            
                           
            self.dicPiList, \
            self.HShannonList, self.HShannonCorrList,\
            self.SeHShannonList, self.SeHShannonCorrList, \
            self.MIlist, self.MIcorrList, self.SeMIList, self.SeMICorrList, self.ijList = \
                self.ent.calcVerticalMutualInfo(desk, self.ent.mySequence.seqs, showmessage=showmessage)
        
            self.entFile.write_VMI_files(sufix, self.dicPiList, self.HShannonList, self.HShannonCorrList,  \
                                         self.SeHShannonList, self.SeHShannonCorrList, \
                                         self.MIlist, self.MIcorrList, \
                                         self.SeMIList, self.SeMICorrList, self.ijList, showmessage=False)
        else:
            stri = "reading MI: %i lines, %i cols, %s"%(len(self.ent.mySequence.seqs), len(self.ent.mySequence.seqs[0]), filename)
            desk.showmsg_obs(stri)

            
            if which_module == 'all':
                ret, msg, self.HShannonList, self.HShannonCorrList = \
                    self.entFile.read_VMI_files(desk.rootTable, sufix, 'hShannon', showmessage=False)
        
                if ret:
                    ret, msg, self.SeHShannonList, self.SeHShannonCorrList = \
                        self.entFile.read_VMI_files(desk.rootTable, sufix, 'seh', showmessage=False)
                        
                ret, msg, self.MIlist, self.MIcorrList = \
                    self.entFile.read_VMI_files(desk.rootTable, sufix, 'mi', showmessage=False)
        
                if ret:
                    ret, msg, self.SeMIList, self.SeMICorrList = \
                        self.entFile.read_VMI_files(desk.rootTable, sufix, 'se', showmessage=False)                        
                
            elif which_module == 'Entropy':
                ret, msg, self.HShannonList, self.HShannonCorrList = \
                    self.entFile.read_VMI_files(desk.rootTable, sufix, 'hShannon', showmessage=False)
        
                if ret:
                    ret, msg, self.SeHShannonList, self.SeHShannonCorrList = \
                        self.entFile.read_VMI_files(desk.rootTable, sufix, 'seh', showmessage=False)
                    
            else:
                ret, msg, self.MIlist, self.MIcorrList = \
                    self.entFile.read_VMI_files(desk.rootTable, sufix, 'mi', showmessage=False)
        
                if ret:
                    ret, msg, self.SeMIList, self.SeMICorrList = \
                        self.entFile.read_VMI_files(desk.rootTable, sufix, 'se', showmessage=False)


            if ret:
                ret, msg, self.dicPiList, self.ijList = \
                    self.entFile.read_VMI_files(desk.rootTable, sufix, 'pij', showmessage=False)


                           
            if not ret:
                print msg + '. Recalculating ....',
                
                self.dicPiList, \
                self.HShannonList, self.HShannonCorrList,\
                self.SeHShannonList, self.SeHShannonCorrList, \
                self.MIlist, self.MIcorrList, self.SeMIList, self.SeMICorrList, self.ijList = \
                    self.ent.calcVerticalMutualInfo(self.desk, self.ent.mySequence.seqs, showmessage=showmessage)

        time.finish()
        print '>>>', time.milli(), 'ms'
        
        return True, 'ok'
    
    def print_Dir_DNA_Histogram(self, title, figNum, listPi, seqZ, lenSeq):
        plt.clf()
        fig = plt.figure()
        # plt.subplot(1, 3, figNum)
        # plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)
        ax = fig.gca(projection='3d')
        
        # cc = lambda arg: colorConverter.to_rgba(arg, alpha=0.6)
        # verts = []
        c = ['r', 'b', 'g', 'y']
        
        for i in range(lenSeq):
            seqY = []
            
            pIs = listPi[i]
            # print pIs
            
            for j in range(4):
                if (pIs[j] > .9) or (pIs[j] < .01):
                    seqY.append(0)
                else:
                    seqY.append(pIs[j])
    
            #print "-----------------"
            #print seqY      
            #print "-----------------"
            # verts.append(list(zip( seqX, seqY) ))
        
            ax.bar(seqZ, seqY, zs=i, zdir='y', color=c, alpha=0.8)
        '''
        poly = PolyCollection(verts, facecolors = [cc('r'), cc('g'), cc('b'), cc('y')])
        poly.set_alpha(0.7)
        ax.add_collection3d(poly, zs=seqZ, zdir='y')
        '''
        ax.set_xlabel('1=A  2=T  3=G   4=C')
        ax.set_xlim3d(0, 4)
        ax.set_ylabel('j (posi seq)')
        ax.set_ylim3d(0, lenSeq)
        ax.set_zlabel('p')
        ax.set_zlim3d(0, 1)
        plt.title("Dirichlet distribution for '" +title+"'")
        
        if figNum == self.maxFig-1:
            plt.show()     
                    
            
    ''' Flavio - updated 07/11/2013 '''       
    def plotHeatMap(self, desk, is3D, MIlist, ijList, L, maxiPos, title, species, limSup, roundVal=4, str_correction=''):
        kL = L-desk.numOfLetters+1
        
        delta = (limSup*1.01)/5.

        if is3D:
            seqX, seqY, seqZ = [[],[],[],[],[]],[[],[],[],[],[]],[[],[],[],[],[]]
            
            for pos in range(len(MIlist)):
                if MIlist[pos] > 0:
                    i,j = ijList[pos]
                    
                    w = int(math.ceil(MIlist[pos]/delta)) - 1
                    if w > 4:
                        w = 4
                        print  'error w > 4', w, delta, limSup
                    try:
                        seqX[w].append(i)
                        seqY[w].append(j)
                        seqZ[w].append(MIlist[pos])
                    except:
                        print  'error w exception', w, delta, limSup
                                       
        else:
            seq = np.zeros( (kL, kL), dtype = float )  # np.zeros((5,), dtype=numpy.int)
            
            ''' i, j, mutual information '''
            for pos in range(len(MIlist)):
                if MIlist[pos] > 0:
                    i,j = ijList[pos]
                    
                    seq[i,j] = MIlist[pos]
                    seq[j,i] = MIlist[pos]
            
            seq = np.array(seq)

        print 'Heat Map valMax', limSup
        tickWidth = 50 if kL < 601 else 100
        ghm = graphPac.HeatMap(desk, is3D, lenSeq=kL, tickWidth=tickWidth)

        stri = ' max(MI) = %xxxf(%xxxf) %s at (x,y)=(%i, %i)'
        stri = stri.replace('xxx','1.'+str(roundVal))
        title2 = stri%(maxiPos[2], maxiPos[3], desk.unit, maxiPos[0], maxiPos[1])
        
        if is3D:
            ghm.plotMI_3D(title+title2, seqX, seqY, seqZ, kL, limSup, desk.unit)
        else:
            ghm.plotMI(title+title2, seq, limSup)

        return ghm


    ''' Flavio - updated 07/11/2013 '''       
    def calcHeatMapDiff(self, listMI_1, listMI_2, ijList, L, title, numOfLetters=1, dna_prot='DNA', limInf=None, limSup=None, showmessage=False):
        ''' rever '''
        
        kL = int(np.ceil(L / numOfLetters))            
        print 'kL',kL
        seq = np.zeros( (kL, kL), dtype = float )  # np.zeros((5,), dtype=numpy.int)
        
        ''' i, j, mutual information '''
        i = 0
        valMax = 0
        valMin = 0
        
        for mat1 in listMI_1:
            mat2 = listMI_2[i]
            i += 1
            diff = mat1[2] - mat2[2]
            seq[mat1[0], mat1[1]] = diff 
            seq[mat1[1], mat1[0]] = diff
            if diff > valMax:
                valMax = diff
            if diff < valMin:
                valMin = diff

             
        ghm = graphPac.HeatMap(lenSeq=len(seq), dna_prot=dna_prot, tickWidth=50)
        

        print valMax, valMin

        ghm.defineMI(title=title, seq=seq,limInf=limInf, limSup=limSup)
    
        return ghm


        
    def print_MI_species_data(self, dicDist, roundVal=4):

        stri = '\n---------------------- JSD(MI) -----------------------------------'
        stri += '\nspecie       \tmean \tSE  \t Conf.Interval(2*SE)'
        stri += '\n-----------------------------------------------------------------------------'
        
        lista = dicDist.keys()
        lista.sort()
        
        for key in lista:
            specie = key
            mi = round(dicDist[key][0],roundVal)
            distError = round(dicDist[key][1],roundVal)

            ci = '[%4.4f, %4.4f]'%(mi-2*distError, mi+2*distError)
            stri += '\n%s \t %4.4f  \t %4.4f  \t %s'%(specie, mi, distError, ci)
            
        return stri + '\n\n'
    

    def build_dna_aa(self, desk):
        if desk.isProtein:
            desk.dna_prot = 'PROT'
            desk.title="Protein Sequence"

            desk.aa = "ACDEFGHIKLMNPQRSTVWY"
            desk.aaList = []
            
            for let1 in desk.aa:
                desk.aaList.append(let1)

            if desk.numOfLetters > 1:
                i = desk.numOfLetters
                
                while (i > 1):
                    i -= 1
                    auxList = desk.aaList
                    desk.aaList = []
    
                    for let1 in auxList:
                        for let2 in desk.aa:
                            desk.aaList.append(let1+let2)
           
            desk.nuc_aa_List = desk.aaList        
            desk.nuc_aa_List = desk.aaList        
        else:
            desk.dna_prot = 'DNA'
            
            desk.protType = 'PROT_ALIGN_CUT'
            desk.title="DNA sequence"

            desk.nuc = "AGTC"
            desk.nucList = []
            
            for let1 in desk.nuc:
                desk.nucList.append(let1)

            if desk.numOfLetters > 1:
                i = desk.numOfLetters
                
                while (i > 1):
                    i -= 1
                    auxList = desk.nucList
                    desk.nucList = []
    
                    for let1 in auxList:
                        for let2 in desk.nuc:
                            desk.nucList.append(let1+let2)
                            
            desk.nuc_aa_List = desk.nucList
            desk.nuc_aa_List = desk.nucList            