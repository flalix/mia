'''
Created on 18/11/2012
Update  on 16/06/2014 -- seq: must pass as string
Update  on 09/10/2014 - MI_Result changed with Standard Errors and MIkCorr
Update  on 09/06/2015 - MI cannot be negative never
Update  on 16/09/2016 - prepareHorizontalMutualInfo_Out

@author: Flavio Lichtenstein
@local: Unifesp Bioinformatica
'''
import BioPythonClass as bioClass
import numpy as np
import classes.Timing as crono
import os

class Mutual_Information_Horizontal:
    def __init__(self, desk, want_nuc_aa_list=False, want_shannon = False):
        self.ok = False

        if want_nuc_aa_list:
            self.build_dna_aa(desk)

        self.pij =  []
        self.MIk = 0
        
        self.time = crono.Timing()
        
        self.entFile = bioClass.Entropy_File()

        self.legColumns = 1
        self.legendTitle = '%s Species'%(desk.organism)
        self.title = 'Horizontal MI'        
        
        ''' numLines=1, numCols=1, 
                       legColumns=self.legColumns, 
                       legendTitle=self.legendTitle, title=self.title '''
        self.ent = bioClass.Entropy(desk)


        
        if want_shannon:
            filenameShannon = 'shannon_random_DNA_Letter%i_Exp100_dic.txt'
            self.filenameShannon1 = filenameShannon%(desk.numOfLetters)
            self.filenameShannon2 = filenameShannon%(desk.numOfLetters*2)
            
            ok, msg = self.ent.read_shannon_random_file(desk.rootEntropy, self.filenameShannon1, desk.numOfLetters, showKeys = False)
            if not ok:
                desk.showmsg_obs(msg + ', didnt find Shannon dictionary for %i letter.'%(desk.numOfLetters))
                return
       
            ok, msg = self.ent.read_shannon_random_file(desk.rootEntropy, self.filenameShannon2, 2*desk.numOfLetters, showKeys = False)
            if not ok:
                desk.showmsg_obs(msg + ', didnt find Shannon dictionary for %i letters.'%(2*desk.numOfLetters))
                return
    
        self.ok = True
        
    def prepareHorizontalMutualInfo(self, desk, sufix, consensus, showmessage=False, hyp_test=False):
        self.desk = desk

        ''' setNames define completeFilename +++ SPECIES !'''
        if not self.entFile.setPrefixSufix(root=desk.rootTable, prefix='HMI', sufix=sufix):
            return

        filename = self.entFile.MI_prefix_sufix_txt(self.entFile.prefix_sufix, 'mi')
        filename = desk.rootTable + filename
        
        if not hyp_test:        
            # print "searching", filename
            self.time.start()
                        
        if not os.path.exists(filename) or desk.de_novo or hyp_test:
            if not hyp_test:
                stri = "running HMI: lines %i cols %i, %s"%(len(self.ent.mySequence.seqs), len(self.ent.mySequence.seqs[0]), filename)
                desk.showmsg_obs(stri)

            self.arrayMI, self.arrayMIcorr, self.arraySE, self.arraySEcorr, self.arrayN =\
                self.run_hmi(desk, self.ent.mySequence.seqs, consensus, hyp_test=hyp_test)
            
            if not hyp_test:
                self.entFile.write_HMI_files(sufix, self.arrayMI, self.arrayMIcorr, self.arraySE, self.arraySEcorr, self.arrayN)
        else:
            stri = "reading HMI: lines %i cols %i, %s"%(len(self.ent.mySequence.seqs), len(self.ent.mySequence.seqs[0]), filename)
            desk.showmsg_obs(stri)
                        
            ok, msg, self.arrayMI, self.arrayMIcorr, _ = \
                self.entFile.read_HMI_files(desk.rootTable, sufix, 'mi', showKeys=False)
            
            if not ok:
                desk.showmsg_obs(msg)
                return False
        
            ok, msg, self.arraySE, self.arraySEcorr, self.arrayN = \
                self.entFile.read_HMI_files(desk.rootTable, sufix, 'se', showKeys=False)
            
            if not ok:
                desk.showmsg_obs(msg)
                return False
        
        if not hyp_test:
            self.time.finish()
            print '>>>', self.time.milli(), 'ms'
        
        return True, 'ok'
 
    ''' resample Ha, get rid from sampledSeqs in H0 '''
    def prepareHorizontalMutualInfo_OneSample(self, desk, sampledSeqs, sufix, consensus, LH0a, showmessage=False):
        self.desk = desk
        ''' setNames define completeFilename +++ SPECIES !'''
        if not self.entFile.setPrefixSufix(root=desk.rootTable, prefix='HMI', sufix=sufix):
            return

        filename = self.entFile.MI_prefix_sufix_txt(self.entFile.prefix_sufix, 'mi')
        filename = desk.rootTable + filename

        self.time.start()
                        
        stri = "run hmi: lines %i cols %i"%(len(self.ent.mySequence.seqs), len(self.ent.mySequence.seqs[0]))
        desk.showmsg_obs(stri)

        self.arrayMI, self.arrayMIcorr, self.arraySE, self.arraySEcorr, self.arrayN =\
            self.run_hmi(desk, self.ent.mySequence.seqs, consensus, LH0a, sampledSeqs)
            
        ''' self.entFile.write_HMI_files(sufix, self.arrayMI, self.arrayMIcorr, self.arraySE, self.arraySEcorr, self.arrayN) '''

        self.time.finish()
        print '>>>', self.time.milli(), 'ms'
        
        return True, 'ok'
 

    def run_hmi(self, desk, seqs, consensus, sampledSeqs=None, hyp_test=False):
        # self.time.start()
        if consensus:
            for i in range(len(seqs)):
                for a in seqs[i]:
                    if a not in ['A','C','G','T']:
                        print "problems with sequence %i letter %s seq="%(i, a, seqs[i])
                        exit()
                        
            self.seqLen  = len(seqs[0])
            self.numSeqs = len(seqs)      
                        
        else:
            ''' if not consensus min(L) can be wrong '''
            minLen = float("inf")
            
            for i in range(len(seqs)):
                stri = ""
                for a in seqs[i]:
                    if a in ['A','C','G','T']:
                        stri += a
                seqs[i] = stri
                # print stri
                
                if len(stri) < minLen:
                    minLen = len(stri)
                    
            # print minLen
            
            self.seqLen = minLen
            self.numSeqs = len(seqs) 
        '''
        seqs = []
        seqs.append("AATTTGGTGATCCTGGATCGTATTGACAATCCTGCTGTCATTGCTGAACTGAAGGCAATCAATCCCAAAGTGACAGTTACCTTCTATCCCTATGATGTTACCGTACCCTTGGCCGAGACAACCAAATTGTTAAACACAATTTTTGCTCAACTTAAAACTGTTGATGTACTGATCAATGGAGCTGGTATTCTTGATGATCATCAGATTGAGCGTACCATTGCAGTTAACTTCACTGGTCTGGTCAATACCACCACCGCCATTTTGGACTTCTGGGATAAGCGTAAGGGTGGTCCAGGTGGTATCATCTGTAACATTGGTTCAGTTACTGGTTTCAATGCCATTTACCAAGTGCCAGTCTATTCTGCATCCAAGGCAGCTGTTGTTAGCTTCACTCAATCCATTGCC")
        seqs.append("AATTTGGTGATCCTGGATCGTATTGACAATGGTGCTGTCATTGCTGAACTGAAGGCAATCAATCCCAAAGTGACAGTTACCTTCTATCCCTATGATGTTACCGCACCCTTGGCCGAGACAACCAAATTGTTAAAGACAATTTTTGCTCAACTTAAGACTGTTGATGTTCTGATGAATGGAGGTGGTATTCTTGATGATCATCAGAATGAGCGTACCATTGGACTGAACTTCACTGGTCTGGTCAATACCACCACCGCCATTTTGGACTTCTGGGATAAGCGTAAGGGTGGTCCAGGTGGTGTCATCTGTAACATTGGTTCAGTTACTGGTTTCAATGCCATTTACCAAGTGCCAGTCTATTCTGCATCCAAGGCAGCTGTTGTTAGCTTCACTCAATCCATTGCC")
        seqs.append("AATTTGGTGATCCTGGATCGTATTGACAATGGTGCTGTCATTGCTGAACTGAGGGCAATCAATCCCAAAGTGACAGTTACCTTCTATCCCTATGATGTTACCGTACCCTTGGCCGAGACAACCAAATTGTTAAAGACAATTTTTGCTCAGCTTAAGACTGTCGATGTTCTGATCAATGGAGCTGGTATTCTTGATGATCATCAGATTGAGCGTACCATTGCAGTTAACTTCACTGGTCTGGTCAATACCACCACCGCCATTTTGGACTTCTGGGATAAGCGTAAGGGTGGTCCAGGTGGTGTCATCTGTAACATTGGTTCAGTTACTGGTTTCAATGCCATTTACCAAGTGCCAGTCTATTCTGCATCCAAGGCAGCTGTTGTTAGCTTCACTCAATCCATTGCC")
        '''
        # self.time.finish()
        # print self.time.milli(), "ms"

        self.calc_Pi(desk, seqs)
        
        arrayMI, arrayMiInf, arraySE, arraySEcorr, arrayN = self.calc_Autocorrelation(desk, seqs, hyp_test)
                
        return arrayMI, arrayMiInf, arraySE, arraySEcorr, arrayN
            
    def calc_Pi(self, desk, seqs):
        self.dicNi = {}
       
        for seq in seqs:
            ''' if aligned seq may vary '''
            len_seq = len(seq)-(desk.offset+desk.numOfLetters-1)
            
            for i in range(desk.offset, len_seq):
                try:
                    self.dicNi[seq[i:i+desk.numOfLetters]] += 1
                except:
                    self.dicNi[seq[i:i+desk.numOfLetters]] = 1
    
        tot = np.sum(self.dicNi.values())
        
        self.dicPi = {}
        for key in self.dicNi.keys():
            self.dicPi[key] = self.dicNi[key] / float(tot)
            
    ''' pass only one seq '''
    def calc_Autocorrelation(self, desk, seqs, hyp_test):
        MIlist, MIcorrList, SElist, SEcorrList, Nlist = [],[],[],[], []
        
        if hyp_test:
            max_k = 151
        else:
            max_k = self.seqLen/2

        ''' from 3 up to L/2 not included '''
        if desk.frame == 0:
            listaK = range(3, max_k, 1)
        elif desk.frame == 1:
            listaK = range(3, max_k, 3)
        elif desk.frame == 2:
            listaK = range(4, max_k, 3)
        else:
            listaK = range(5, max_k, 3)

        for k in listaK:
            if not hyp_test:
                if k%50 == 0:
                    desk.showmsg_obs(k, same_line = True)
            # print '-->k',k
            # MIk, pij, pi, qj = self.autocorrelation_K_without_correction(desk, seq, k, self.nucList, offset, numOfLetters)
            ''' correction: 09-13/10/2014 '''
            MIk, MIkinf, SE, SEcorr, N = self.autocorrelation_K(desk, seqs, k)

            if MIk < 0:
                desk.showmsg_obs('>>> MI is negative in calc_Autocorrelation for k=%i  MI = %f"'%(k, MIk) )
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
        
        if not hyp_test:
            desk.showmsg_obs('.')
        return MIlist, MIcorrList, SElist, SEcorrList, Nlist 

              
    def autocorrelation_K(self, desk, seqs, k):
        '''
            frame 0: k = 3,4,5,6,7,8,9,10 ...
            frame 1: k = 3,6,9 ...
            frame 2: k = 4,7,10 ...
            frame 3: k = 5,8,11 ...
        '''
        dicPij = {}
        dicCij = {}
        
        
        ''' self.nuc_aa_List == list Nuc or AA '''

        N = 0
        for seq in seqs:
            ''' if aligned len can vary '''
            lenTot = len(seq)-k-desk.numOfLetters-1
            
            ''' SUM OF 3 FRAMES '''
            for ini in range(0, lenTot): 
                nij = seq[ini:ini+desk.numOfLetters]+ seq[ini+k:ini+k+desk.numOfLetters]
                if len(nij) != desk.numOfLetters*2:
                    pass
                N += 1
                
                try:
                    dicCij[nij] += 1
                except:
                    dicCij[nij] = 1
        
        if not dicCij:
            pass
        
        tot = float(np.array(dicCij.values()).sum())
        
        for nij in dicCij.keys():
            dicPij[nij] = dicCij[nij] / tot

        ''' dicCi -- marginal frequency of i in ij '''
        dicPi = {}
        dicCi = {}
        
        for key1 in self.dicNi.keys():
            tot = 0
            for key2 in self.dicNi.keys():
                try:
                    tot += dicCij[key1+key2]
                except:
                    pass
                
            dicCi[key1] = tot

        tot = float( np.sum(dicCi.values()) )
        for key1 in dicCi.keys():
            try:
                dicPi[key1] = dicCi[key1] / tot
            except:
                pass
                
        ''' how many i states exist '''
        Bi = sum([1 for val in dicCi.values() if val > 0])
        
        ''' dicCj -- marginal frequency of j in ij '''
        dicQj = {}
        dicCj = {}
        for key2 in self.dicNi.keys():
            tot = 0
            for key1 in self.dicNi.keys():
                try:
                    tot += dicCij[key1+key2]
                except:
                    pass
                
            dicCj[key2] = tot
            
        tot = float( np.sum(dicCj.values()))
        if tot==0:
            print "Error tot==0, dicCj=", dicCj
            print "dicNi=", self.dicNi
            exit()
            
        for key in dicCj.keys():
            dicQj[key] = dicCj[key] / tot
                
        ''' how many j states exist '''
        Bj = sum([1 for val in dicCj.values() if val > 0])
        
        ''' how many ij states exist '''
        Bij = sum([1 for val in dicCij.values() if val > 0])

        ''' Roulston bias correction '''
        correction = (Bi+Bj-Bij-1) / (2.*N)
        # print 'Bij %i,  Bi %i  Bj %i, and N = %i, correction=%f'%(Bij, Bi, Bj, N, correction)

        totMI = 0.
        totVar = 0
        totVarCorr = 0

        for key in dicPij.keys():
            n1 = key[0:desk.numOfLetters]
            n2 = key[desk.numOfLetters:desk.numOfLetters+desk.numOfLetters]
    
            Pij = dicPij[key]

            if Pij > 0 and Pij < 1:
                Pi = dicPi[n1]
                Qj = dicQj[n2]
                mi = Pij * np.log(Pij/(Pi*Qj))

                totVarCorr +=  (np.log(Pi) + np.log(Qj) - np.log(Pij) + mi)**2 * (Pij*(1.-Pij))
                totVar += mi**2 * (Pij*(1.-Pij))
                totMI += mi

        if totMI < 0:            
            print("*** Impossible ***totMI = %f, for k =%i"%(totMI, k) )
            raw_input('Enter to continue.')

        
        '''   MI cannot be negative never: 09/06/2015  '''
        totMIcorr = totMI + correction
        if totMIcorr < 0:
            # print('Bias correction turn to negative !', totMIcorr, totMI, correction)
            totMIcorr = 0

        
        return totMI, totMIcorr, np.sqrt(totVar/N), np.sqrt(totVarCorr/N), N

    def build_dna_aa(self, desk):
        ''' from IUPAC data 
        protein_letters = "ACDEFGHIKLMNPQRSTVWY"
        extended_protein_letters = "ACDEFGHIKLMNPQRSTVWYBXZJUO"
        unambiguous_dna_letters = "GATC"
        '''
        
        if desk.isProtein:
            # http://biopython.org/DIST/docs/api/Bio.PDB.Polypeptide-pysrc.html
            self.aa = "ACDEFGHIKLMNPQRSTVWY"
            self.aaList = []
            
            for let1 in self.aa:
                self.aaList.append(let1)

            if desk.numOfLetters > 1:
                i = desk.numOfLetters
                
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

            if desk.numOfLetters > 1:
                i = desk.numOfLetters
                
                while (i > 1):
                    i -= 1
                    auxList = self.nucList
                    self.nucList = []
    
                    for let1 in auxList:
                        for let2 in self.nuc:
                            self.nucList.append(let1+let2)

            self.nuc_aa_List = self.nucList

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
   
   
    