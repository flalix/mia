'''
Created on 24/11/2014

@author: Dabreu
'''
import numpy as np
import classes.BioPythonClass as bioClass
import classes.Timing as crono
import classes.BarGraphic as graphPac
import os
import matplotlib.pyplot as plt
from scipy import stats

class test_Class:
    def __init__(self, isProtein=False, numOfLetters=1):
        self.isProtein = isProtein
        self.numOfLetters = numOfLetters
        self.offset = 0

        ''' from IUPAC data 
        protein_letters = "ACDEFGHIKLMNPQRSTVWY"
        extended_protein_letters = "ACDEFGHIKLMNPQRSTVWYBXZJUO"
        unambiguous_dna_letters = "GATC"
        
        http://biopython.org/DIST/docs/api/Bio.PDB.Polypeptide-pysrc.html
        '''
        
        if self.isProtein:
            self.aa = "ACDEFGHIKLMNPQRSTVWY"
            self.aaList = []
            
            for let1 in self.aa:
                self.aaList.append(let1)

            if self.numOfLetters > 1:
                i = self.numOfLetters
                
                while (i > 1):
                    i -= 1
                    auxList = self.aaList
                    self.aaList = []
    
                    for let1 in auxList:
                        for let2 in self.aa:
                            self.aaList.append(let1+let2)
           
            self.nuc_aa_List = self.aaList
            
        else:
            self.nuc = "AGTC"
            self.nucList = []
            
            for let1 in self.nuc:
                self.nucList.append(let1)

            if self.numOfLetters > 1:
                i = self.numOfLetters
                
                while (i > 1):
                    i -= 1
                    auxList = self.nucList
                    self.nucList = []
    
                    for let1 in auxList:
                        for let2 in self.nuc:
                            self.nucList.append(let1+let2)
                            
            self.nuc_aa_List = self.nucList
                    
    def run_hmi(self, seqs):
        self.seqLen  = len(seqs[0])
        self.numSeqs = len(seqs)
        
        self.calc_Pi(seqs)
        
        arrayMi, arrayMiInf, arraySE, arraySEcorr, arrayN = self.calc_Autocorrelation(seqs)
        
        return arrayMi, arrayMiInf, arraySE, arraySEcorr, arrayN
            
    def calc_Pi(self, seqs):
        self.dicNi = {}
        
        len_seq = len(seqs[0])-(self.offset+self.numOfLetters-1)
        
        for seq in seqs:
            for i in range(self.offset, len_seq):
                try:
                    self.dicNi[seq[i:i+self.numOfLetters]] += 1
                except:
                    self.dicNi[seq[i:i+self.numOfLetters]] = 1
    
        tot = np.sum(self.dicNi.values())
        
        self.dicPi = {}
        print '--- calc Pi'
        for key in self.dicNi.keys():
            self.dicPi[key] = self.dicNi[key] / float(tot)
            print '%s = %f'%(key,self.dicPi[key])        
        print '--------------'
        
    ''' pass only one seq '''
    def calc_Autocorrelation(self, seqs):
        MIlist, MIcorrList, SElist, SEcorrList, Nlist = [],[],[],[], []
        
        max_k = self.seqLen/2
        print 'max_k = %i'%(max_k)

        if self.frame == 0:
            listaK = range(3, max_k, 1)
        elif self.frame == 1:
            listaK = range(3, max_k, 3)
        elif self.frame == 2:
            listaK = range(4, max_k, 3)
        else:
            listaK = range(5, max_k, 3)

        print 'listaK=', listaK
        
        for k in listaK:
            if k%50 == 0:
                print k
            # print '-->k',k
            # MIk, pij, pi, qj = self.autocorrelation_K_without_correction(desk, seq, k, self.nucList, offset, numOfLetters)
            ''' correction: 09-13/10/2014 '''
            MIk, MIkinf, SE, SEcorr, N = self.autocorrelation_K(seqs, k)

            if MIk < 0:
                print '>>> MI is negative in calc_Autocorrelation for k=%i  MI = %f"'%(k, MIk)
                MIk = 0
                
            '''
                MI = mutual information (horizonta)
                k = k distance for HMI
                MIkinf = MIk infinite = bias correction
                varMI = Var[MI]
                N = total elements searched
            '''
            MIlist.append(MIk)
            MIcorrList.append(MIkinf)
            SElist.append(SE)
            SEcorrList.append(SEcorr)
            Nlist.append(N)
                        
        return MIlist, MIcorrList, SElist, SEcorrList, Nlist 
 
 
    def autocorrelation_K(self, seqs, k):
        '''
            frame 0: k = 3,4,5,6,7,8,9,10 ...
            frame 1: k = 3,6,9 ...
            frame 2: k = 4,7,10 ...
            frame 3: k = 5,8,11 ...
        '''
        dicPij = {}
        dicCij = {}
        lenTot = self.seqLen-k-(self.numOfLetters -1)
        
        ''' self.nuc_aa_List == list Nuc or AA '''

        N = 0
        print '\n--- k=%i --------------------------------'%(k)
        
        for seq in seqs:
            ''' SUM OF 3 FRAMES '''
            print '--->',seq,": ",
            for ini in range(0, lenTot): 
                nij = seq[ini:ini+self.numOfLetters]+ seq[ini+k:ini+k+self.numOfLetters]
                print nij,
                N += 1
                
                try:
                    dicCij[nij] += 1
                except:
                    dicCij[nij] = 1
                    
                #print dicCij[nij],
            print ''        
        
        tot = float(np.array(dicCij.values()).sum())
        
        for nij in dicCij.keys():
            dicPij[nij] = dicCij[nij] / tot
            print 'kmer=%s, count=%i, perc=%f'%(nij, dicCij[nij], dicPij[nij])
            
        print '---- tot = %f, sum=1 == %i'%(float(tot), sum(dicPij.values()))

        ''' dicCi -- marginal frequency of i in ij '''
        dicPi = {}
        dicCi = {}
        print '--- first kmer -----'
        for key1 in self.dicNi.keys():
            tot = 0
            for key2 in self.dicNi.keys():
                try:
                    tot += dicCij[key1+key2]
                except:
                    continue
                
            dicCi[key1] = tot
            print "'%s' = %i"%(key1, dicCi[key1])


        tot = float( np.sum(dicCi.values()) )
        for key1 in dicCi.keys():
            dicPi[key1] = dicCi[key1] / tot
                
        ''' how many i states exist '''
        Bi = sum([1 for val in dicCi.values() if val > 0])
        print 'len marginal i=', Bi, 'dicCi.values()', dicCi.values()
        
        ''' dicCj -- marginal frequency of j in ij '''
        dicQj = {}
        dicCj = {}
        print '--- second kmer -----'
        for key2 in self.dicNi.keys():
            tot = 0
            for key1 in self.dicNi.keys():
                try:
                    tot += dicCij[key1+key2]
                except:
                    pass
                
            dicCj[key2] = tot
            print key2, dicCj[key2]
            
        tot = float( np.sum(dicCj.values()))
        for key in dicCj.keys():
            dicQj[key] = dicCj[key] / tot
        
        ''' how many j states exist '''
        Bj = sum([1 for val in dicCj.values() if val > 0])
        
        print 'len marginal j=', Bj, 'dicCj.values()', dicCj.values()
        
        ''' how many ij states exist '''
        Bij = sum([1 for val in dicCij.values() if val > 0])

        ''' Roulston bias correction '''
        correction = (Bi+Bj-Bij-1) / (2.*N)
        print 'Bij %i,  Bi %i  Bj %i, and N = %i, correction=%f'%(Bij, Bi, Bj, N, correction)

        totMI = 0.
        totVar = 0
        totVarCorr = 0
        
        print 'nuc1,\tnuc2,\tPij,\t\tPi,\t\tQj,\t\tmi,\t\tMIk'
        for key in dicPij.keys():
            n1 = key[0]
            n2 = key[1]
    
            Pij = dicPij[key]

            if Pij > 0 and Pij < 1:
                Pi = dicPi[n1]
                Qj = dicQj[n2]
                mi = Pij * np.log(Pij/(Pi*Qj))

                totVarCorr +=  (np.log(Pi) + np.log(Qj) - np.log(Pij) + mi)**2 * (Pij*(1.-Pij))
                totVar += mi**2 * (Pij*(1.-Pij))
                totMI += mi
                print '%s \t%s \t%f \t%f \t%f \t%f \t%f'%(n1, n2, Pij, Pi, Qj, mi, totMI)

                    
        if totMI < 0:
            print("*** Impossible ***MIk = %f, for k =%i"%(totMI, k) )
            raw_input('Enter to continue.')

        
        ''' Groeesse et al. 2000 - MIk = 1/3 * sum(Pi*Pj) all 3 frames
        MIk = MIk/3.
        varMI = np.sqrt(totVar / float(N))
        varMICorr = np.sqrt(totVarCorr / float(N))  '''


        return totMI, totMI + correction, np.sqrt(totVar/N), np.sqrt(totVarCorr/N), N


class Vertical_MI:
    def __init__(self, desk, maxFig=1, showmessage=False):
        self.basic = bioClass.Basic()
        self.ok = False
        self.desk = desk
        self.names = desk.names
        
        if not self.names:
            desk.showmsg_obs('Please give the names of the sequences')
            return
        
        self.left = 0.05
        self.top = 0.75
        
        self.isProtein = desk.isProtein

        self.organism   = desk.organism
        self.gene       = desk.gene
        self.species    = desk.species
        self.seqType    = desk.seqType

        self.cutoffLength = desk.cutoffLength
        self.cutoffNumSeq = desk.cutoffNumSeq
        self.numOfLetters = desk.numOfLetters

        self.root = desk.rootFasta
        self.rootImage = desk.rootImage
        self.rootTable = desk.rootTable
        
        filenameShannon = 'shannon_random_DNA_Letter%i_Exp100_dic.txt'
        self.filenameShannon1 = filenameShannon%(1)
        self.filenameShannon2 = filenameShannon%(2)

        self.names = desk.names
        self.isProtein = desk.isProtein
        
        
        if desk.isProtein:
            self.dna_prot = 'PROT'
            self.title="Protein Sequence"
        else:
            self.dna_prot = 'DNA'
            
            self.protType = 'PROT_ALIGN_CUT'
            self.title="DNA sequence"

        
        self.numSpecies = len(self.names)
        self.manySpecies = []
        
        self.maxFig = maxFig

        self.species = []
        self.entFile = bioClass.Entropy_File()

        self.legColumns = 1
        self.legendTitle = 'Drosophila Species'

        self.ent = bioClass.Entropy(root=self.root, filename='dummy', read=False, numOfLetters=self.numOfLetters,
                       numLines=1, numCols=1, 
                       left=self.left, top=self.top, legColumns=self.legColumns, 
                       legendTitle=self.legendTitle, title=self.title, showmessage=showmessage, 
                       alfa= None, instanceGraph=False)


        ok, msg = self.ent.read_shannon_random_file(self.rootTable, self.filenameShannon1, 1, showKeys = False)
        if not ok:
            desk.showmsg_obs(msg + ', didnt find Shannon dictionary for 1 letter.')
            return
   
        ok, msg = self.ent.read_shannon_random_file(self.rootTable, self.filenameShannon2, 2, showKeys = False)
        if not ok:
            desk.showmsg_obs(msg + ', didnt find Shannon dictionary for 2 letters.')
            return
    
        ''' logp = log(D).mean(axis=0) eq 05 - 
            testar bem o param para separar as curvas '''
        param = 0.001
        hum = 1.00000
        zero = 0.00000
        
        self.tA = np.array(self.noZeros([hum,zero,zero,zero], param))
        self.tT = np.array(self.noZeros([zero,hum,zero,zero], param)) 
        self.tG = np.array(self.noZeros([zero,zero,hum,zero], param)) 
        self.tC = np.array(self.noZeros([zero,zero,zero,hum], param)) 
        '''

        self.tA = np.array([1.,0.,0.,0.])
        self.tT = np.array([0.,1.,0.,0.])  
        self.tG = np.array([0.,0.,1.,0.])
        self.tC = np.array([0.,0.,0.,1.])       
        '''
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
        
    def read_file(self, tName, showmessage=False):
        
        self.ent.mySequence.filename = tName[0]
        self.ent.mySequence.name = tName[1]

        print self.ent.mySequence.filename 
        print 'name', self.ent.mySequence.name
        
        # load in sequence
        root_filename_fasta = self.root + self.ent.mySequence.filename
        
        try:
            if not self.ent.mySequence.read_fasta(root_filename_fasta):
                print 'could not read fasta filename:' + self.ent.mySequence.filename
                return False
        except:
            print 'Error occured, reading fasta filename:' + self.ent.mySequence.filename
            return False
        
        self.manySpecies.append(self.ent.mySequence.sequence)
        return True
        
    def read_file_new(self, tName, showmessage=False):
        
        self.ent.mySequence.filename = tName[0]
        self.ent.mySequence.name = tName[1]

        if not self.ent.mySequence.read_fasta_new(self.root + self.ent.mySequence.filename):
            print 'could not read %s'%(self.root + self.ent.mySequence.filename)
            return False
        
        return True
        

    def prepareVerticalMutualInfo(self, desk, iLoop, sufix, showmessage=False):
        self.desk = desk
        numOfLetters = desk.numOfLetters
        isProtein = desk.isProtein
        de_novo = desk.de_novo
        
        '''
        try:
            print 'first seq', self.ent.mySequence.sequence[0][:50]
        except:
            print 'first seq', self.ent.mySequence.sequence[0]
        '''    
        print 'Prepare Mutual Information: num indiv: %i length: %i' % \
              (self.ent.mySequence.getNumOfSeqs(),self.ent.mySequence.getLengthSequence())
        
        ''' setNames define completeFilename +++ SPECIE !'''
        if not self.entFile.setPrefixSufix(root=self.rootTable, prefix='VMI', sufix=sufix):
            return

        filename = self.entFile.MI_prefix_sufix_txt(self.entFile.prefix_sufix, 'mij')
        filename = self.rootTable + filename

        tim = crono.Timing()
        tim.start()

        if not os.path.exists(filename) or de_novo:
            stri = "new: lines %i cols %i, %s"%(len(self.ent.mySequence.seqs), len(self.ent.mySequence.seqs[0]), filename)
            print stri
            
            ''' MIlist, MIcorrList, listIJ, SElist '''
            self.MIlist, self.MIcorrList, self.listIJ, self.SElist = \
                self.ent.calcVerticalMutualInfo(self.desk, self.ent.mySequence.seqs, showmessage=showmessage)
        
            
            # if not os.path.exists(filename):
            self.entFile.write_VMI_files(sufix, self.MIlist, self.MIcorrList, self.listIJ, self.SElist, showmessage=False)
        else:
            ''' first MI file '''
            stri = "reading MI: %i lines, %i cols, %s"%(len(self.ent.mySequence.seqs), len(self.ent.mySequence.seqs[0]), filename)
            print stri
            
            ret, msg, self.MIlist, self.MIcorrList = \
                self.entFile.read_VMI_files(self.rootTable, sufix, 'mij', showmessage=False)

            if not ret:
                return ret, msg
            
            ret, msg, self.SElist, self.listIJ = \
                self.entFile.read_VMI_files(self.rootTable, sufix, 'var', showmessage=False)

            if not ret:
                return ret, msg
            
            if not ret:
                print msg + '\nRecalculating ....'
                
                self.MIlist, self.MIcorrList, self.listIJ, self.SElist = \
                    self.ent.calcVerticalMutualInfo(self.desk, self.ent.mySequence.seqs, numOfLetters=numOfLetters, isProtein=isProtein, showmessage=showmessage)

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
    def calcHeatMap(self, MIlist, unit, listIJ, L, maxiPos, title2, gene, cutoffLength, cutoffNumSeq, seqType, numOfLetters=1, dna_prot='DNA', showmessage=False,specie='', limInf=None, limSup=None, roundVal=4, str_correction=''):
        title = title2[0]
        merged = title2[1]
        
        kL = int(np.ceil(L / numOfLetters))            
        print 'kL',kL
        seq = np.zeros( (kL, kL), dtype = float )  # np.zeros((5,), dtype=numpy.int)

        valMax = 0
        valMin = 0  # float('Inf')
               
        ''' i, j, mutual information '''
        for pos in range(len(MIlist)):
            i,j = listIJ[pos]
            
            seq[i,j] = MIlist[pos]
            seq[j,i] = MIlist[pos]

            if MIlist[pos] > valMax:
                valMax = MIlist[pos]
            if MIlist[pos] < valMin:
                valMin = MIlist[pos]

        if not limSup:
            limSup = valMax
        elif valMax > limSup:
            limSup = valMax

        if not limInf:
            limInf = valMin
        elif valMin < limInf:
            limInf = valMin

        limInf = 0
        
        print 'Heat Map valMax', valMax
        
        gmi = graphPac.HeatMap(lenSeq=len(seq), dna_prot=dna_prot, tickWidth=50)
        
        if specie == '':
            title = title + ' ' +  self.organism
        else:
            title = title + ' ' + self.organism + ' ' + specie
           
        try: 
            stri = 'max(MI) = %xxxf(%xxxf) %s at (x,y)=(%i, %i)'
            stri = stri.replace('xxx','1.'+str(roundVal))
            title2 = stri%(maxiPos[2], maxiPos[3], unit, maxiPos[0], maxiPos[1])
        except:
            print 'except ', maxiPos
            title2 = ''
            
        title2 += '\nGene %s, %s seqs %s, cutoffLength %i, cutoffNumSeq %i #letter %i' % (gene, merged, seqType, cutoffLength, cutoffNumSeq, numOfLetters) 

        gmi.plotMI(title=title, seq=seq, numOfSeqs=self.ent.mySequence.getNumOfSeqs(),limInf=limInf, limSup=limSup, title2=title2)
    
        return gmi

    ''' Flavio - updated 07/11/2013 '''       
    def calcHeatMapDiff(self, MIlist_1, MIlist_2, listIJ, L, title, numOfLetters=1, dna_prot='DNA', limInf=None, limSup=None, showmessage=False):
        ''' rever '''
        
        kL = int(np.ceil(L / numOfLetters))            
        print 'kL',kL
        seq = np.zeros( (kL, kL), dtype = float )  # np.zeros((5,), dtype=numpy.int)
        
        ''' i, j, mutual information '''
        i = 0
        valMax = 0
        valMin = 0
        
        for mat1 in MIlist_1:
            mat2 = MIlist_2[i]
            i += 1
            diff = mat1[2] - mat2[2]
            seq[mat1[0], mat1[1]] = diff 
            seq[mat1[1], mat1[0]] = diff
            if diff > valMax:
                valMax = diff
            if diff < valMin:
                valMin = diff

             
        gmi = graphPac.HeatMap(lenSeq=len(seq), dna_prot=dna_prot, tickWidth=50)
        

        print valMax, valMin

        gmi.defineMI(title=title, seq=seq,limInf=limInf, limSup=limSup)
    
        return gmi

    
    def write_file(self, root, filename, stri):
        filename = root + filename
        try:
            f = open(filename, 'w')
            s = str(stri)
            f.write(s)
            f.flush()
            
            print 'write file %s'%filename
            
        except:
            print 'Could not write %s'%filename
            return False
        
        finally:
            f.close()
        
        
        return True
    
    def show_JSD_MI(self, dicDist, factor, title, xlabel, ylabel = 'JSD', showGraph=True, saveGraph=False, roundVal=4, root="C:/data/image/", filename="my_graph", fileType='png', isLog=False, sigmas=1):
       
        title = title + ' (error=%i*SE)'%(sigmas)
        i = 0
        lista = dicDist.keys()
        lista.sort()
        
        seqX = []
        seqY = []
        seqError = []
        seqSeq = []

        has_negative = False
        for key in lista:
            if isLog:
                try:
                    d = np.array(dicDist[key][0])
                    se = np.array(dicDist[key][1])
                    
                    dist = np.round(np.log10(d),roundVal)
               
                    if not has_negative and dist < 0:
                        has_negative = True
                except:
                    dist = 0


                try:
                    distError = np.abs( np.round(np.log10(d)     ,roundVal) - \
                                        np.round(np.log10(d+se)) ,roundVal)
                    
                except:
                    distError = 0
            else:
                d = np.array(dicDist[key][0])
                se = np.array(dicDist[key][1])
                
                if factor == 1:
                    dist      = np.round(d,roundVal)
                    distError = np.round(se,roundVal)
                else:
                    dist      = np.round(d*factor,roundVal)
                    distError = np.round(se*factor,roundVal)
                
            distError *= sigmas
            
            seqX.append(i)
            seqSeq.append(key)
            seqY.append(dist)
            seqError.append(distError)
            i += 1

        
        '''

        _, ax = plt.subplots()
        ax.bar(seqX, seqY, yerr=seqError,color='y')
        ax.set_xticklabels( seqSeq )
        ax.set_xticks(seqX)
        '''
        plt.close("all")
        plt.clf()            
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        ax.bar(seqX, seqY, yerr=seqError,color='y')
        ax.set_xticklabels(seqSeq)
        ax.set_xticks(seqX)

            
        if has_negative:
            plt.ylim(min(seqY), max(seqY))
        

        plt.xticks(rotation=90)
        plt.subplots_adjust(bottom=.2) 
        plt.tick_params(axis='both', which='major', labelsize=8)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)



        if showGraph:
            plt.show()
        
        if saveGraph:
            filename = root + filename + '.' + fileType
            plt.savefig(filename, format=fileType)
        
        
    def calc_anova(self, seqMI):
            
        #one way anova
        f_value, p_value = stats.f_oneway(*seqMI)
        
        print '\n-----------'
        if p_value <= 0.05:
            stri = 'At least one distribution is statically different from the others, by ANOVA'
        else:
            stri =  'The distributions are statiscally similar, by ANOVA'
            
        stri += '\nf_value %f   p_value %.5e' %(f_value, p_value)
        
        return stri + '\n\n'
        
    def print_MI_species_data(self, dicDist, roundVal=4):

        stri = '\n---------------------- JSD(MI) -----------------------------------'
        stri += '\nspecie       \tmean \terror  \t Conf.Interval'
        stri += '\n-----------------------------------------------------------------------------'
        
        lista = dicDist.keys()
        lista.sort()
        
        for key in lista:
            specie = key
            mi = round(dicDist[key][0],roundVal)
            distError = round(dicDist[key][1],roundVal)

            ci = '[%4.4f, %4.4f]'%(mi-distError, mi+distError)
            stri += '\n%s \t %4.4f  \t %4.4f  \t %s'%(specie, mi, distError, ci)
            
        return stri + '\n\n'
    
    def print_data_params(self, speciesParams, roundVal=4, filename=None, stri = '', saveData=False):
        str_param = stri + '\n'
        str_param += 'species \t#seqs \tmean \tSdv \tmedian \tConf.Interval 95% \tmaximum \tminimum\n'

        for mat in speciesParams:
            ''' seqSpecies.append([species,sp,numOfSeqs, meanVal, medianVal, stdVal]) '''
            species = mat[0]
            numOfSeqs = mat[2]
            mean = round(mat[3],roundVal)
            median = round(mat[4],roundVal)
            sdv = round(mat[5],roundVal)
            maxi = round(mat[6],roundVal)
            mini = round(mat[7],roundVal)
            infLim = round(mean - 1.96 * sdv / np.sqrt(numOfSeqs),roundVal)
            supLim = round(mean + 1.96 * sdv / np.sqrt(numOfSeqs),roundVal)
            
            stri = '%s \t%i \t%.4f \t%.4f \t%.4f \t[%.4f, %.4f] \t%.4f \t%.4f\n'
            stri = stri.replace('%.4f', '%.'+str(roundVal)+'f')
            str_param += stri%(species, numOfSeqs, mean, sdv, median, infLim, supLim, maxi, mini) 


        # print str_param
        if saveData:
            self.write_file(self.rootTable, filename, str_param)
