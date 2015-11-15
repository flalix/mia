#-*- coding: utf-8 -*-
'''
Created on 18/06/2012
Updated on 25/06/2012
Updated on 08/05/2013  read_fasta
Updated on 12/02/2014  countManyLettersDic, calcShannonVerticalEntropy
Update  on 09/06/2015 - MI cannot be negative never
Update  on 14/09/2015 - read_fasta_samples 
Update  on 27/09/2015 - full randomizing and shuffling: read_fasta_samples(..., sampledSeqs=None, method=None)

@author: Flavio Lichtenstein
@local: Unifesp Bioinformatica
'''
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
# from scipy.stats import f_oneway 
from scipy.stats import f_oneway
from operator import itemgetter #, attrgetter
import copy, os, numpy as np
from math import pow
import random
# from numpy.random.mtrand import shuffle
     
class MySequence:
    def __init__(self, desk, showmessage=False):
        self.clearSequences()
        
        self.NucleotideAlphabet = ["A","C","G","T"]

        
    def read(self, desk, filename, showmessage):
        self.name = filename.replace('.fasta','')
        self.name = filename.replace('.fas','')
        self.filename = filename

        self.root_filename_fasta = self.root + self.filename
        
        self.fileNameAln = desk.root_fasta + self.name + '.aln'
        self.fileNameMscAln = desk.root_fasta + self.name + '_Muscle.aln'
        self.fileNameDnd = desk.root_fasta + self.name + '.dnd'
        self.fileOutputAln = desk.root_fasta + self.name + ' - aligned.fasta'
        self.fileOutputAlnCut = desk.root_fasta + self.name + ' and cut.fasta'

        self.fileOutputName = self.filename + ' - aligned and cut.fasta' 
        self.fileOutputProteinName = self.filename + ' - protein - aligned and cut.fasta'
        
        self.numSequences = 0
        
        self.basic = Basic()
        
        self.seqs = []
        self.seqNCBI = self.readMySequence(desk.isProtein, showmessage)

        if self.seqNCBI:
            for j in range(len(self.seqNCBI)):
                self.seqs.append(self.seqNCBI[j].seq)
                # print 'self.seqNCBI[j].seq', self.seqNCBI[j].seq
                        
            self.numSequences = len(self.seqs)
            self.lenSequence = len(self.seqs[0])
        else:
            self.numSequences = 0
            self.lenSequence = 0           

        
    def clearSequences(self):
        self.seqs = []
        self.seqCut = []
        self.seqFinal = []
        self.idList = []
        self.seqNCBI = []
        self.seqProteinAligned = []
            
    def set_seq_all_param(self, desk, species, seqs, which, setSample=None):
        if seqs:
            self.seqs = seqs
            
            self.numSequences = len(self.seqs)
            self.numIndiv = self.numSequences
            self.lenSequence = len(self.seqs[0])
        else:
            self.seqs = []
            
            self.numSequences = 0
            self.numIndiv = 0
            self.lenSequence = 0
            
        str_sample = "_sample" if setSample else ""
        
        desk.which = which
        
        if desk.which == 'VMI' or desk.which == 'VSH':
            self.sufix_vmi = ('%s_%s_%s_%s_%s_NOL%i_%iL_cutoff%i%s%s') %\
                            (desk.organism, desk.minmax, species, desk.seqType, desk.gene, desk.numOfLetters, \
                             desk.cutoffLength, desk.cutoffNumSeq, str_sample, desk.stri_random)
        else:
            ''' hmi_Drosophila_americana_Gene_frame0_NOL1_400L_cutoff10_MI.txt '''
            self.sufix_hmi = ('%s_%s_%s_%s_%s_frame%i_NOL%i_%iL_cutoff%i%s%s') %\
                             (desk.organism, desk.minmax, species, desk.seqType, desk.gene, desk.frame, \
                              desk.numOfLetters, desk.cutoffLength, desk.cutoffNumSeq, str_sample, desk.stri_random)


    def set_seq(self, seqs):
        self.seqs = seqs
        # print '**** ', self.seqs
        self.numSequences = len(self.seqs)
        self.lenSequence = len(self.seqs[0])
        
        
        
    ''' https://github.com/biopython/biopython/blob/master/Doc/examples/fasta_dictionary.py '''
    def get_accession_num(self, seq_record):
        self.accession_atoms = seq_record.id.split('|')
        self.gb_name = self.accession_atoms[3]
        # strip the version info before returning
        return self.gb_name[:-2]
        
    def readMySequence(self, isProtein, showmessage):
        seq = []
        
        if isProtein:
            proteinCompl = ' - protein.fasta'
            self.completeName = self.filename.replace('.fasta','') + proteinCompl
        else:
            self.completeName = self.filename
            
        print "'readMySequence', reading ", self.completeName

        if (os.path.exists(self.completeName)):
            try:
                for seq_record in SeqIO.parse(self.completeName, "fasta"):
                    seq.append(seq_record)
                    
                if (showmessage):
                    print "'" + self.completeName + "' num seq: ", len(seq)
    
            except ValueError:
                print "Exception: error reading file ''" + self.completeName + "''"

        else:
            print "Error reading file ''" + self.completeName + "''"
       
        return seq


    def findInNCBI(self, arrProtId, i):
        
        name = 'dummy'
        print name, 'len', len(arrProtId)-1, 'ids:', arrProtId, '\n'
        seqNCBI = []

        Entrez.email = self.email
        # handle = Entrez.efetch(db="nucleotide", rettype="fasta", id="6273291")

        for j in range(len(arrProtId)):
            print arrProtId[j]
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", id=arrProtId[j])
            
            seqNCBI.append(SeqIO.read(handle, "fasta"))
                              
            handle.close()
            
            print "%s with %i features" % (seqNCBI[j].id, len(seqNCBI[j].features))
            print seqNCBI[j], '\n'

        return seqNCBI
   
        
    def save_fasta_dna_prot(self, seqParam, dna_prot, filename):
        my_records = []
        
        for i in range(len(seqParam)):
            stri_seq = ''
            stri_head = ''

            for j in range(len(seqParam[i].seq)):
                stri_seq += seqParam[i].seq[j]

            for j in range(len(seqParam[i].id)):
                stri_head += seqParam[i].id[j]

            if (dna_prot == 'D'):
                rec = SeqRecord(Seq(stri_seq, IUPAC.unambiguous_dna),
                                id=stri_head, description='')
            else:
                rec = SeqRecord(Seq(stri_seq),
                    id=stri_head, description='')
                
            my_records.append(rec)
    
        # saving the file           
        try:
            output_handle = os.open(filename, "w")
            SeqIO.write(my_records, output_handle, "fasta")
            print "File '" + filename + "' saved."
            ret = True
            
        except ValueError:
            print "File '" + filename + "' not saved. Writing error: ", ValueError
            ret = False

        finally:    
            output_handle.close()

        return ret
    
    def save_fasta_add(self, filename, seqs, which="fasta"):
        if which not in ['clustal','emboss','fasta', 'fasta-m10','ig','maf','nexus','phylip','phylip-sequential', 'phylip-relaxed', 'stockholm']:
            print 'Which is the format?'
            return False
        
        try:
            handle = open(filename, "a")
            SeqIO.write(seqs, handle, which)
            ret = True
            
        except ValueError:
            print "File '" + filename + "' not saved. Writing error: ", ValueError
            ret = False

        finally:    
            handle.close()
            
        return ret

    def save_fasta_new(self, filename, seqs):
        
        try:
            handle = open(filename, "w")
            SeqIO.write(seqs, handle, 'fasta')
            print "Saved '%s'."%(filename)
            ret = True
        except:
            print "File '%s' not saved. Writing error."%(filename)
            ret = False
        finally:
            handle.close()
            
        return ret
    

    def save_fasta(self, seq, filename):
        print '\nsaving protein in fasta:' + filename

        my_records = []
        
        for i in range(self.getNumSequences()):
            stri_aa = ''
            stri_head = ''
            
            striProt = seq[i]
            striHead = self.idList[i]
            
            # print striProt
            
            for j in range(len(striProt)):
                stri_aa += striProt[j]

            for j in range(len(striHead)):
                stri_head += striHead[j]
            
            # print striHead[j], striProt[j]
            rec = SeqRecord(Seq(stri_aa), id=stri_head, description='')
            # print rec
                
            my_records.append(rec)
    

        try:
            if os.path.isfile(filename):
                os.remove(filename)
                print "\n----------------------------------------"
                print "Deleting file '" + filename + "'"
            print "----------------------------------------"
        except ValueError:
            print "File '" + filename + "'could not be deleted", ValueError
            return False
                        
        try:
            # output_handle = os.open(filename, os.O_RDWR|os.O_CREAT )
            output_handle = open(filename, "w")
            print "\n----------------------------------------"
            print "File '" + filename + "'",
            SeqIO.write(my_records, output_handle, "fasta")
            print " saved."
            print "----------------------------------------"

        except ValueError:
            print "File '" + filename + "' not saved. Writing error: ", ValueError
            return False
        
        return True


    def save_fasta2(self, seqGi, seq, seqDesc, filename):
        # filename=filename.replace('.fasta','')
        print '\nsaving protein in fasta2:' + filename

        my_records = []
        
        if not seqDesc:
            seqDesc = ''
        
        self.numSequences = len(seq)
        print 'self.getNumSequences()', self.getNumSequences()

        for i in range(self.getNumSequences()):
            # print seq[i], seqGi[i], seqDesc[i]
            rec = SeqRecord(Seq(seq[i]), id=seqGi[i], description=seqDesc[i])
                
            my_records.append(rec)
    

        print 'filename', filename
        try:
            if os.path.isfile(filename):
                os.remove(filename)
        except ValueError:
            print "File '" + filename + "'could not be deleted", ValueError
            return False
                        
        try:
            # os.open does not works
            output_handle = open(filename, "w")
            # output_handle = os.open(filename, os.O_EXCL | os.O_CREAT)
            SeqIO.write(my_records, output_handle, "fasta")
            print "File '" + filename + "' saved."

        except ValueError:
            print "File '" + filename + "' not saved. Writing error: ", ValueError
            return False
        
        return True
    

    def read_fasta(self, root_filename_fasta, showmessage=False):
        self.root_filename_fasta = root_filename_fasta
       
        self.seqs = []
        
        try:
            SeqIO.parse(self.root_filename_fasta, format="fasta").next()
        except:
            try:
                self.root_filename_fasta = self.root_filename_fasta.replace(".fasta",".fas")
                SeqIO.parse(self.root_filename_fasta, format="fasta").next()
            except:
                self.numSequences = 0
                self.lenSequence = 0
                print "** Error reading fasta file ''" + root_filename_fasta + "''"
                return False
                    
        for seq_record in SeqIO.parse(self.root_filename_fasta, format="fasta"):
            self.add(seq_record.seq, seq_record.id)
            
        self.numSequences = len(self.seqs)
        self.lenSequence = len(self.seqs[0])

        if showmessage:
            print 'Num of indiv=%i, length=%i, read %s'%(self.numSequences, len(self.seqs[0]), root_filename_fasta)

        return True

    ''' Flavio 14/09/2015'''
    def read_fasta_samples(self, root_filename_fasta, method, sampledSeqs=None, showmessage=False):
        self.root_filename_fasta = root_filename_fasta
        
        self.seqs = []
        
        try:
            SeqIO.parse(self.root_filename_fasta, format="fasta").next()
        except:
            try:
                self.root_filename_fasta = self.root_filename_fasta.replace(".fasta",".fas")
                SeqIO.parse(self.root_filename_fasta, format="fasta").next()
            except:
                self.numSequences = 0
                self.lenSequence = 0
                print "** Error reading fasta file ''" + self.root_filename_fasta + "''"
                return False

        if method == "shuffle":
            for seq_record in SeqIO.parse(self.root_filename_fasta, format="fasta"):
                seq_record = self.shuffling(seq_record)
                self.add(seq_record.seq, seq_record.id)
        else:
            ''' randomize: sampling seqs or not '''
            if sampledSeqs:
                i = 0
                ''' randomizing with known records, or pure random '''
                for seq_record in SeqIO.parse(self.root_filename_fasta, format="fasta"):
                    if i in sampledSeqs:
                        self.add(seq_record.seq, seq_record.id)
                        
                    i += 1
            else:
                ''' full randomizing with the same length '''
                for seq_record in SeqIO.parse(self.root_filename_fasta, format="fasta"):
                    L = len(seq_record.seq)
                    seq = [self.NucleotideAlphabet[random.randint(0, 3)] for _ in range(L)]
                    seq = "".join(seq)
                    seq_record = SeqRecord(Seq(seq, IUPAC.unambiguous_dna), id="random", description='random')
                    self.add(seq_record.seq, seq_record.id)
                
                
        self.numSequences = len(self.seqs)
        self.lenSequence = len(self.seqs[0])

        if showmessage:
            print 'Num of indiv=%i, length=%i, read %s'%(self.numSequences, len(self.seqs[0]), root_filename_fasta)

        return True


    def shuffling(self, x):
        if type(x)==str:
            x = [a for a in x]
            random.shuffle(x)   
            x = ''.join(str(e) for e in x)   
            return x
        
        if type(x)==int:
            if x < 0:
                isNegative = True
                x = -x
            else:
                isNegative = False
                
            x = [a for a in str(x)]
            random.shuffle(x)   
            x = ''.join(str(e) for e in x) 
            return -int(x) if isNegative else int(x)
        
        if type(x)==float:
            if x < 0:
                isNegative = True
                x = -x
            else:
                isNegative = False
                
            x = [a for a in str(x)]
            random.shuffle(x)   
            x = ''.join(str(e) for e in x)   
            return -float(x) if isNegative else float(x)  
        
        if type(x)==list:
            random.shuffle(x)   
            return x
    
        if type(x)==Seq:
            seq = str(x)
            seq = [x for x in seq]
            random.shuffle(seq)   
            seq = ''.join(str(e) for e in seq)  
            x = Seq(seq, IUPAC.unambiguous_dna)
            return x
    
        if type(x)==SeqRecord:
            idi = x.id
            desc = x.description
            name = x.name
    
            seq = str(x.seq)
            seq = [x for x in seq]
            random.shuffle(seq)   
            seq = ''.join(str(e) for e in seq)  
            x = Seq(seq, IUPAC.unambiguous_dna)
            return SeqRecord(x, id=idi, description=desc, name=name)
        
        print("Could not shuffle")
        return x
        

    def fasta_file_simulation(self, desk, L, N):
        self.filename = "Simulation"
        self.name = "Simulation"
        self.root_filename_fasta = desk.rootFasta + self.filename
    
        self.seqs = []
        for i in range(N):
            seq = [random.sample("ACGT", 1)[0] for _ in range(L)]
            seq = ''.join(str(e) for e in seq)
            self.seqs.append(seq)

            seq_record = SeqRecord(Seq(seq, IUPAC.unambiguous_dna),
                                id="Simulation"+str(i), description='')
            self.add(seq_record.seq, seq_record.id)

            
        self.numSequences = N
        self.lenSequence = L
        
        print 'Num of indiv=%i, length=%i, read %s'%(self.numSequences, len(self.seqs[0]), self.filename)

        return True

    def read_file_one_sample(self, desk, filename, name, isReference, species, N, sampledSeqs):
        
        self.filename = filename
        self.name = name

        self.root_filename_fasta = desk.rootFasta + self.filename
        
        self.seqs = []
        
        try:
            SeqIO.parse(self.root_filename_fasta, format="fasta").next()
        except:
            try:
                self.root_filename_fasta = self.root_filename_fasta.replace(".fasta",".fas")
                SeqIO.parse(self.root_filename_fasta, format="fasta").next()
            except:
                self.numSequences = 0
                self.lenSequence = 0
                print "** Error reading fasta file ''" + self.root_filename_fasta + "''"
                return False
     
        i = 0
        for seq_record in SeqIO.parse(self.root_filename_fasta, format="fasta"):
            ''' if different species: add all - but we should read HMI ...'''
            if species != desk.refSpecies:
                self.add(seq_record.seq, seq_record.id)
            else:
                if isReference:
                    ''' if reference, sampledSeqs are the sampled'''
                    if i in sampledSeqs: self.add(seq_record.seq, seq_record.id)
                else:
                    ''' not reference, the same species, exclude sampledSeqs '''
                    if i not in sampledSeqs: self.add(seq_record.seq, seq_record.id)
                
                i += 1

            
        self.numSequences = len(self.seqs)
        self.lenSequence = len(self.seqs[0])

        print 'Num of indiv=%i, length=%i, read %s'%(self.numSequences, len(self.seqs[0]), self.root_filename_fasta)

        return True

    def read_file_only_samples(self, desk, filename, name, isReference, sampledSeqs, method):
        
        self.filename = filename
        self.name = name

        self.root_filename_fasta = desk.rootFasta + self.filename
        
        try:
            SeqIO.parse(self.root_filename_fasta, format="fasta").next()
        except:
            try:
                self.root_filename_fasta = self.root_filename_fasta.replace(".fasta",".fas")
                SeqIO.parse(self.root_filename_fasta, format="fasta").next()
            except:
                self.numSequences = 0
                self.lenSequence = 0
                desk.showmsg_obs("** Error reading fasta file ''" + self.root_filename_fasta + "''")
                return False
     
        i = 0
        self.seqs = []

        if not method or method == "standard":
            if isReference:
                ''' this sample '''
                for seq_record in SeqIO.parse(self.root_filename_fasta, format="fasta"):
                    ''' if IS reference, sampledSeqs are the samples'''
                    if i in sampledSeqs: self.add(seq_record.seq, seq_record.id)
                    i += 1
            else:
                ''' complement sample '''
                for seq_record in SeqIO.parse(self.root_filename_fasta, format="fasta"):
                    ''' if NOT reference: complement sampledSeqs'''
                    if i not in sampledSeqs: self.add(seq_record.seq, seq_record.id)
                    i += 1
        elif method == "shuffle":
            if isReference:
                ''' this sample '''
                for seq_record in SeqIO.parse(self.root_filename_fasta, format="fasta"):
                    ''' if IS reference, sampledSeqs are the samples'''
                    if i in sampledSeqs:
                        seq_record = self.shuffling(seq_record)
                        self.add(seq_record.seq, seq_record.id)
            
                    i += 1
            else:
                ''' complement sample '''
                for seq_record in SeqIO.parse(self.root_filename_fasta, format="fasta"):
                    ''' if NOT reference: complement sampledSeqs'''
                    if i not in sampledSeqs: 
                        seq_record = self.shuffling(seq_record)
                        self.add(seq_record.seq, seq_record.id)
            
                    i += 1
        else:
            ''' randomize: sampling seqs or not '''
            if isReference:
                ''' this sample '''
                for seq_record in SeqIO.parse(self.root_filename_fasta, format="fasta"):
                    ''' if IS reference, sampledSeqs are the samples'''
                    if i in sampledSeqs:
                        seq_record = self.shuffling(seq_record)
                        L = len(seq_record.seq)
                        seq = [self.NucleotideAlphabet[random.randint(0, 3)] for _ in range(L)]
                        seq = "".join(seq)
                        seq_record = SeqRecord(Seq(seq, IUPAC.unambiguous_dna), id="random", description='random')
                        self.add(seq_record.seq, seq_record.id)
            
                    i += 1
            else:
                ''' complement sample '''
                for seq_record in SeqIO.parse(self.root_filename_fasta, format="fasta"):
                    ''' if NOT reference: complement sampledSeqs'''
                    if i not in sampledSeqs: 
                        seq_record = self.shuffling(seq_record)
                        L = len(seq_record.seq)
                        seq = [self.NucleotideAlphabet[random.randint(0, 3)] for _ in range(L)]
                        seq = "".join(seq)
                        seq_record = SeqRecord(Seq(seq, IUPAC.unambiguous_dna), id="random", description='random')
                        self.add(seq_record.seq, seq_record.id)
            
                    i += 1
                    
               
        if not self.seqs:
            desk.showmsg_obs('Error: No seqs selected - total loop=%i, proposed=%s - %s'%(i, str(sampledSeqs),filename))
            return False
        
        self.numSequences = len(self.seqs)
        self.lenSequence = len(self.seqs[0])

        ''' RFOS = read file only samples '''
        desk.showmsg_obs('Sampling=%i/%i, length=%i - %s - %s'%(self.numSequences, i, len(self.seqs[0]), str(sampledSeqs),filename))
        

        return True

    def read_fasta_3methods(self, desk, filename, name, method=None, sampledSeqs=None, showmessage=False):
        self.filename = filename
        self.name = name

        root_filename = desk.rootFasta + self.filename

        if not method or method == "" or method == "standard":
            if not self.read_fasta(root_filename, showmessage):
                return False 
        else:
            if not self.read_fasta_samples(root_filename, method, sampledSeqs, showmessage):
                return False
        
        return True

    
        
    def read_fasta_file(self, desk, tName, showmessage=False):
        
        self.ent.mySequence.filename = tName[0]
        self.ent.mySequence.name = tName[1]

        print self.ent.mySequence.filename 
        
        # load in sequence
        root_filename_fasta = desk.rootFasta + self.ent.mySequence.filename
        
        try:
            if not self.ent.mySequence.read_fasta(root_filename_fasta):
                print 'could not read fasta filename:' + self.ent.mySequence.filename
                return False
        except:
            print 'Error occured, reading fasta filename:' + self.ent.mySequence.filename
            return False

        return True
        
    
    def split_species(self, gene, organism,showmessage=False):
        dicSpecies = {}
        
        ''' organism|species|gene|Prod|product 
            195116780_XM_002002894_Drosophila|mojavensis|DmojAdh1|Prod|Adh1-PA
        '''
        
        for numSeq in range(self.getNumSequences()):
            
            seq  = self.seqs[numSeq]
            desc = self.idList[numSeq]
        
            pos = desc.find('_')
            off = 1
            if pos < 0:
                print 'bad structure - description'
                return None, None
                
            # seqId = desc[:pos]
            _ = desc[:pos]
            stri = desc[pos+off:]
                
            pos = stri.find('_')
            off = 1
            if pos < 0:
                print 'bad structure - accession number'
                return None, None

            ''' see http://www.ncbi.nlm.nih.gov/refseq/about/  
                like: XM_002002894 '''
            if pos < 5:
                ncbi_prefix = stri[:pos+off]

                stri = stri[pos+off:]
                ''' find the second _ '''
                pos = stri.find('_')

                if pos < 0:
                    print 'bad structure - accession number'
                    return None, None

                
                #ncbi_prefix = ncbi_prefix + stri[:pos] # accNum
                _ = ncbi_prefix + stri[:pos] # accNum
            else:
                _ = stri[:pos] # accNum
                
            stri = stri[pos+off:]
            # print seqId, accNum
            
            ''' organism|species|gene|Prod|product 
                195116780_XM_002002894_Drosophila|mojavensis|DmojAdh1|Prod|Adh1-PA
            '''

            pos = stri.find('|')
            off = 1
            if pos < 0:
                print 'bad structure - species'
                return None, None

            # organism = stri[:pos]
            _ = stri[:pos]
            stri = stri[pos+off:]
                
            pos = stri.find('|')
            off = 1
            if pos < 0:
                print 'bad structure - species'
                return None, None

            species = stri[:pos]
            species = species.strip()
            species = species.replace(" ","_")
        
            if species not in dicSpecies.keys():
                dicSpecies[species] = []
        
            lista = dicSpecies[species]
            lista.append([seq, desc])

        return dicSpecies, sorted(list(dicSpecies.keys()))        


    def split_species_seqs(self, organism, gene, seqs, showmessage=False):
        dicSpecies = {}
        
        ''' organism|species|gene|Prod|product 
            195116780_XM_002002894_Drosophila|mojavensis|DmojAdh1|Prod|Adh1-PA
        '''
        
        for numSeq in range(len(seqs)):
            
            seq  = seqs[numSeq].seq
            idSeq = seqs[numSeq].id
        
            pos = idSeq.find('_')
            off = 1
            if pos < 0:
                print 'bad structure - description'
                return None, None
                
            stri = idSeq[pos+off:]
                
            pos = stri.find('_')
            off = 1
            if pos < 0:
                print 'bad structure - accession number'
                return None, None

            ''' see http://www.ncbi.nlm.nih.gov/refseq/about/  
                like: XM_002002894 '''
            if pos < 5:
                ncbi_prefix = stri[:pos+off]

                stri = stri[pos+off:]
                ''' find the second _ '''
                pos = stri.find('_')

                if pos < 0:
                    print 'bad structure - accession number'
                    return None, None

                
                #ncbi_prefix = ncbi_prefix + stri[:pos] # accNum
                _ = ncbi_prefix + stri[:pos] # accNum
            else:
                _ = stri[:pos] # accNum
                
            stri = stri[pos+off:]
            # print seqId, accNum
            
            ''' organism|species|gene|Prod|product 
                195116780_XM_002002894_Drosophila|mojavensis|DmojAdh1|Prod|Adh1-PA
            '''

            pos = stri.find('|')
            off = 1
            if pos < 0:
                print 'bad structure - species'
                return None, None

            # organism = stri[:pos]
            _ = stri[:pos]
            stri = stri[pos+off:]
                
            pos = stri.find('|')
            off = 1
            if pos < 0:
                print 'bad structure - species'
                return None, None

            species = stri[:pos]
            species = species.strip()
            species = species.replace(" ","_")
        
            if species not in dicSpecies.keys():
                dicSpecies[species] = []
        
            lista = dicSpecies[species]
            lista.append([seq, idSeq])

        return dicSpecies, sorted(list(dicSpecies.keys())) 
    
    def strangeAAchar(self, seq, showmessage):
        if (showmessage):
            print '\nsearching for AA strange characters'
        
        problem = False
        
        for j in range(len(seq)):
            
            stri = seq[j]
            
            for i in range(len(stri)):
                choose = stri[i]
                
                if (Basic.getStringAA().find(choose) < 0):
                    problem = True
                    if (showmessage):
                        print 'Strange char "' + choose + '" seq ', j, ' pos', i
                    else:
                        break
    
        return problem 
        
            
    def strangeNucelotideChar(self, seq, showmessage):
        if (showmessage):
            print '\nsearching for strange characters'
        
        problem = False
        stringNuc = self.basic.getDnaNucleotideString()
        
        for j in range(len(seq)):
            stri = seq[j]
            
            for i in range(len(stri)):
                choose = stri[i]
                
                if (stringNuc.find(choose) < 0):
                    problem = True
                    if (showmessage):
                        print 'Strange char "' + choose + '" seq ', j, ' pos', i
                    else:
                        break
    
        return problem        
    
    def add(self, seq, _id):
        self.seqs.append(seq)
        self.idList.append(_id)

        
    def addProtein(self, protein):
        self.seqs.append(protein)

    def getNumSequences(self):
        return self.numSequences
            
    def getLengthSequence(self):
        if (self.seqs == None):
            return 0
        return len(self.seqs[0])

    def getLengthSeqCut(self):
        if (self.seqCut == None):
            return 0
        return len(self.seqCut[0])
        
    def getLengthSeqFinal(self):
        return len(self.seqFinal)
        
    def getNumOfSeqs(self):
        if (self.seqs == None):
            return 0
        
        return len(self.seqs)
      
    def getNumOfProteinSeqs(self):
        if (self.seqs == None):
            return 0
        
        return len(self.seqs)
            
    def getSequence(self, i):
        return self.seqs[i]
    
    def getSequenceCut(self, i):
        return self.seqCut[i]
    
    def getSequenceFinal(self, i):
        return self.seqFinal[i]
    
    def getProtein(self, i):
        return self.seqs[i]
    
    def writeCutSequence(self, i, stri):
        self.seqCut[i] = stri
    
    def writeFinalSequence(self, i, stri):
        self.seqFinal[i] = stri
    
    def cutBeginSeq(self, seq, posIni=-1, showmessage=False):
        if posIni > 0:
            posLeft = posIni-1
        else:
            posLeft = 0
            
            for j in range(len(seq)):
                
                stri = seq[j]
    
                for i in range(len(stri)):
                    if (stri[i] == '-'): 
                        pos = i
                    else:
                        break
                        
                    if (pos != 0):
                        if (pos > posLeft):
                            posLeft = pos

        print posLeft
        for j in range(len(seq)):
            
            stri = seq[j]
            
            self.seqCut.append(stri[posLeft + 1:])
            
            if (showmessage):
                print 'left cut ', self.seqCut[j]    
            
        return posLeft+1, self.seqCut
             
                 
    def cutEndSeq(self, seq, posEnd=-1, showmessage=False):
        if posEnd > 0:
            posRight = posEnd+1
        else:
            posRight = 99999999999
            
            for j in range(len(seq)):
                stri = seq[j]
        
                i = len(stri) - 1;        
                while True:
                    # print 'begin lenght ', i
                    if (stri[i] == '-'): 
                        pos = i
                    else:
                        break
                        
                    if (pos != 0):
                        if (pos < posRight):
                            posRight = pos
                            
                    i = i - 1
                    if (i < 0):
                        break;

        print posRight
        for j in range(len(seq)):

            stri = seq[j][:posRight]
            self.writeCutSequence(j, stri)
            
            if (showmessage):
                print 'rigth cut ', self.seqCut[j]
                
        return self.seqCut

    
    def cutAAStrangeChar(self, stringAA, showmessage):
        posLeft = -1
        
        if (len(self.seqs) == 0):
            return
            
        found = 0
        for j in range(self.getnumSequences()):
            stri = self.seqs[j]
            
            for i in range(len(stri)):
                choose = stri[i] 
                if (stringAA.find(choose) < 0):
                    print stringAA, 'char', choose
                    found = j
                    posLeft = i
                    break
                
            if (posLeft >= 0):
                break
                        
        if (posLeft == -1):
            return False
        
        if (showmessage):
            print 'cutting strange chars ', len(self.seqs[0]), ' found at ', found, 'char', choose
            print 'x' * (posLeft + 4), '  --> pos ', posLeft, '  length: ', len(self.seqs[found]), 'char', choose
            if (j < 2): 
                print 'ini', self.seqs[found]

        for j in range(self.getnumSequences()):
            stri = self.seqs[j]

            stri = stri[:posLeft] + stri[posLeft + 1:]
            if (showmessage):
                if (j < 2): 
                    print 'end', stri
                    
            self.writeProtein(j, stri)

    
        return True
    
    def cutAminoAcid(self, position, showmessage):
        
        for j in range(self.getnumSequences()):
            stri = self.seqs[j]
            stri = stri[:position] + stri[position + 1:]
            
            if (showmessage):
                if (j < 5): 
                    print stri

            self.writeProtein(j, stri)
    
        if (showmessage):
            print 'x' * position, '  --> pos ', position, '  lenght: ', len(self.seqs[0])
    
        return True
    
    def cutMiddleSeq(self, seqCut, isProtein=True, showmessage=False):
        posLeft = -1
        
        if isProtein:
            stringVocab = self.basic.getStringAA()
        else:
            stringVocab = self.basic.getDnaNucleotideString()
            
        found = 0
        for j in range(len(seqCut)):
            stri = seqCut[j]

            for i in range(len(stri)):
                choose = stri[i]
                if (stringVocab.find(choose) < 0):
                    posLeft = i
                    found = j
                    break
                
            if (posLeft >= 0):
                break                
                        
        if (posLeft == -1):
            return False, seqCut
        
        if (posLeft % 3 == 0):
            offsetLeft = 0
            offsetRight = 3
            
        elif (posLeft % 3 == 1):
            offsetLeft = 1
            offsetRight = 2
        else:
            offsetLeft = 2
            offsetRight = 1
                 
        if (showmessage):
            print 'cutting middle width ', len(seqCut[0]), ' found at ', found, 'char', choose
            print 'x' * (posLeft + 4), '  --> pos ', posLeft, ' length: ', len(seqCut[found]), 'char', choose
            if (j < 2): 
                print 'ini', seqCut[found]

        for j in range(self.getNumSequences()):
            stri = seqCut[j]

            stri = stri[:posLeft - offsetLeft] + stri[posLeft + offsetRight:]
            if (showmessage):
                if (j < 2): 
                    print 'end', stri

            seqCut[j] = stri

    
        return True, seqCut

    def convertProtein(self):
        
        print 'convertProtein'
        for j in range(self.numSequences):
            
            dna = ''
            stri = self.seqFinal[j]
            
            for i in range(len(stri)):
                dna += stri[i]
    
            try:
                coding_dna = Seq(dna, IUPAC.unambiguous_dna)
            except ValueError:
                print "Error Verifying Seq: ", ValueError
                return False
   
            try:
                protein = coding_dna.translate()
            except ValueError:
                print "Error Translanting: ", ValueError
                return False
                        
            # print '  ',j,'. self: ',protein 
            self.addProtein(protein)
            
        return True


    def transpose(self, seqs, showmessage):
        setTransp = []
        
        if (showmessage):
            print '----------- Transposing --------------'
            print 'width ', len(seqs[0]), 'num seqs ', len(seqs)
        
        # print 'transposing ... seqs ', len(self.seqs)
        # j = columns
        for col in range(len(seqs[0])):
            stri = ''
            
            # i = lines
            for row in range(len(seqs)):
                striRow = seqs[row]
                # print '--> ' + striRow
                stri += striRow[col]
            
            setTransp.append(stri)
            '''
            if (showmessage):
                print "column: ", col, ". ",setTransp[col]
            '''
        if (showmessage):
            print '----------- Transposing End --------------'
            print 'width ', len(setTransp[0]), 'num seq ', len(setTransp)
            
        return setTransp 
    
    def build_DNA_consensus(self, seqs, listPolymorphism, force=False):
        '''     ['A','G','T','C','?'] '''
        for posi in listPolymorphism:
            result = np.array([0,0,0,0,0])
            
            for seq in seqs:
                if seq[posi] == 'A':
                    result += np.array([1,0,0,0,0])
                elif seq[posi] == 'G':
                    result += np.array([0,1,0,0,0])
                elif seq[posi] == 'T':
                    result += np.array([0,0,1,0,0])
                elif seq[posi] == 'C':
                    result += np.array([0,0,0,1,0])
                else:
                    result += np.array([0,0,0,0,1])
            # print posi, result
            
            dna = ['A','G','T','C']
            possible = []
            # soma = sum(result)
            
            for i in range(len(result)-1):
                possible.append( [dna[i], result[i]] )

            
            possible = sorted(possible, key=itemgetter(1), reverse=True)
       
            # print 'possible', possible
            converted = True

            for i in range(len(seqs)):
                if seqs[i][posi] not in dna:
                    NUC = seqs[i][posi]

                    for maxi in possible:
                        converted = True
                        
                        # print "  .pos:%i trying to convert '%s' having %s in %s"%(posi, NUC, str(maxi), str(possible)), 
                        if NUC == 'R':
                            #print 'R = Purine (A or G)',
                            if maxi[0] not in ['A','G']:
                                converted = False
                                
                        elif NUC == 'Y':
                            #print 'Y = Pyrimidine (C, T, or U)',
                            if maxi[0] not in ['C','T','U']:
                                converted = False
                                
                        elif NUC == 'M':
                            #print 'M = C or A',
                            if maxi[0] not in ['A','C']:
                                converted = False

                        elif NUC == 'K':
                            #print 'K = T, U, or G',
                            if maxi[0] not in ['G','T','U']:
                                converted = False

                        elif NUC == 'W':
                            #print 'W = T, U, or A',' x ', maxi[0]
                            if maxi[0] not in ['A','T','U']:
                                converted = False

                        elif NUC == 'S':
                            #print 'S = C or G',
                            if maxi[0] not in ['C','G']:
                                converted = False

                        elif NUC == 'B':
                            #print 'B = C, T, U, or G (not A)',
                            if maxi[0] not in ['C','T','U','G']:
                                converted = False

                        elif NUC == 'D':
                            #print 'D = A, T, U, or G (not C)',
                            if maxi[0] not in ['A','T','U','G']:
                                converted = False

                        elif NUC == 'H':
                            #print 'H = A, T, U, or C (not G)',
                            if maxi[0] not in ['A','T','U','C']:
                                converted = False

                        elif NUC == 'V':
                            #print 'V = A, C, or G (not T, not U)',
                            if maxi[0] not in ['A','C','G']:
                                converted = False

                        elif NUC == 'N':
                            #print 'N = any base (A, C, G, T, or U)',
                            if maxi[0] not in ['A','C','G','T','U']:
                                converted = False
                            
                        elif NUC == '-':
                            converted = True
                            
                        else:
                            print NUC, 'What is this ?????'
                            converted = False

                        if converted:
                            # print " >> changing '%s' with '%s'"%(NUC, maxi[0])
                            seqs[i]  = seqs[i][0:posi] + maxi[0] + seqs[i][posi+1:len(seqs[i])]
                            break
                        #else:
                            #print 'not converted in', maxi[0]

                    if not converted:
                        if force:
                            seqs[i]  = seqs[i][0:posi] + possible[0][0] + seqs[i][posi+1:len(seqs[i])]
                            print ' >> Forcing ', NUC, 'with', possible[0][0]
                        else:
                            print ' !error! did not change ', NUC, ' in', possible


        return seqs
        
class Basic(object):
    def __init__(self):
        self.log202 = np.log10(20) / np.log10(2)
        self.log2 = np.log10(2)
        self.ln2 = np.log(2)
        self.inv_ln2 = 1./self.ln2
        self.ln2_2 = self.ln2 ** 2
        self.convShannonTsallis = 1.417199516  # =1/LN(2)
    
    def getDnaNucleotides(self):
        return ['A','T','G','C']
    
    def getDnaNucleotideString(self):
        return 'ATGC'
    
    def getRnaNucleotides(self):
        return ['A','U','G','C']
    
    def getrnaNucleotideString(self):
        return 'AUGC'

    # in Hydropathy index[95] http://en.wikipedia.org/wiki/Amino_acid
    def getSeqAA(self):
        return ['A', 'M', 'C', 'F', 'L', 'V', 'I', 'G', 'T', 'S', 'W', 'Y', 'P', 'H', 'N', 'D', 'E', 'Q', 'K', 'R']
    
    def getStringAA(self):
        return 'AMCFLVIGTSWYPHNDEQKR'
    
    def getAaPos(self, string, aa):
        return string.find(aa)
    

    def sum(self, seq):
        x = 0.
        
        for i in range(len(seq)):
            x += seq[i]
        
        return x
    

    def find_root(self, q, Sq):
        p = .5
        Smin = (1. - (pow(p,q) + (pow((1-p),q) )) ) / (q - 1.) 
        
        if Smin > Sq:
            print 'error  Smin > Sq,   Smin=', Smin, '  Sq', Sq
            return -1 

        p = .499
        pb = 0  # pbefore
        found = False
        
        # print '\n --------------------------------'
        # print '-----  q=', q, ' Sq=', Sq
        # print ' --------------------------------'
        
        Scalc = (1. - (pow(p,q) + (pow((1-p),q) )) ) / (q - 1.)  
        
        dif = Scalc - Sq
        
        error = 0.01
        
        delP = 0.025
        change_signal = False
        i = 0
        
        while True:
            if change_signal:
                delP = abs(p - pb)*.43
                # print '    *** change signal '

            if dif < -error:
                isNegative = True
                pant = p
                p -= delP
                
                if p < 0.000000001:
                    pant = delP
                    p = delP / 2.3
                    delP = p
                    
                pb = pant
                # print '    --> negative  '
        
            elif dif > error:
                isNegative = False
                pant = p
                p += delP
                
                if p >= .5:
                    delP = 2*delP
                    p = pant - delP
                    
                    if p <= 0:
                        delP = delP / 4.3
                        p = delP
                        pant = 2. * delP
                    
                # print '    --> positive  '
                pb = pant

        
            else:
                # print '         ___found p =', p, 'dif', dif
                found = True
                break
            
            i += 1
            if i > 100:
                # print 'count off', 'p', p, 'dif', dif, 'i', i
                return p
        
            if (p == 0):
                Scalc = 0
            else:
                Scalc = (1. - (pow(p,q) + (pow((1-p),q) )) ) / (q - 1.)  

            dif = Scalc - Sq
            # print i, 'p', p, 'dif', dif, ' Scalc', Scalc, 'dif=', dif, ' q', q

            if dif > 0:
                if isNegative:
                    change_signal = True
                else:
                    change_signal = False
            else:
                if isNegative:
                    change_signal = False
                else:
                    change_signal = True 
                    
                      
        if not found:
            # print '******* Not Found'
            return 1.
        
        # print '******* root =', p
        return p
    
    
    def tira_aspas(self, stri):
        if (stri == None):
            return 'Null'
        
        try:
            stri = stri.replace("'","")
            stri = stri.replace('"',"'")
        except:
            stri = str(stri)
            stri = stri.replace("'","")
            stri = stri.replace('"',"'")
            stri = stri.replace('[',"")
            stri = stri.replace(']',"")
            
        return stri 
    
class Entropy_File(object):
    def __init__(self):
        self.root = ''
        self.filename = ''
        self.completeFilename = ''
        
        self.shannon_dic = None
        self.tsallis_dic = None


    def setPrefixSufix(self, root, prefix, sufix):
        self.root = root
        self.prefix = prefix
        self.sufix = sufix
        self.prefix_sufix = prefix + '_' + sufix
        self.complete_Prefix_sufix = root + self.prefix_sufix
        self.fileExists = False
        
        if (not root) or (root == ''):
            print 'Define root and filename'
            return False
                    
        if not os.path.exists(root):
            os.mkdir(root)
            
            if not os.chdir(root):
                print 'could not create root' + root
                return False
            
        return True

    
    def setNames(self, root, filename):
        self.root = root
        self.filename = filename
        self.completeFilename = root + filename
        self.fileExists = False
        
        if (not root) or (root == '') or (not filename) or (filename ==''):
            print 'Define root and filename'
            return False
                    
        if not os.path.exists(root):
            os.mkdir(root)
            
            if not os.chdir(root):
                print 'could not create root' + root
                return False
            
        self.fileExists = os.path.isfile(self.completeFilename)
        
        return True
        


    def read_HMI_files(self, root, sufix, which, showKeys = False):
        filename = root + 'HMI_%s_%s.txt'%(sufix, which)
        if which == "mi":
            print "\tReading HMI", filename
       
        if not os.path.isfile(filename):
            print 'read_HMI_files error', filename
            return False, "Could not find: '" + filename +"'", None, None, None
                    
        try:
            f = open(filename, 'r')
            horiz_mi_dic = eval(f.read())
            f.close()
            
        except IOError, e:
            print 'error', e  # e.errno, e.strerror
            return False, "Could not open: '" + self.completeFilename +"'", None, None, None
    
        if showKeys:
            print '\n-------------------------------' 
            print '  reading keys' 
            keys = horiz_mi_dic.keys()
            
            print 'keys =', keys, ' - ', type(keys)
            print '-------------------------------' 

        if which == 'mi':
                return True, "HMI MIlist Table Opened " + filename, \
                    horiz_mi_dic['MI'], horiz_mi_dic['MIcorr'], None

        if which == 'se':
                return True, "HMI SE Table Opened " + filename, \
                    horiz_mi_dic['SE'], horiz_mi_dic['SEcorr'], horiz_mi_dic['N']

        return True, "HMi pij Table Opened " + filename, \
                    horiz_mi_dic['Pij'], horiz_mi_dic['Pi'], None


    def read_VMI_files(self, root, sufix, which, showmessage=False):
        ret, msg, listA, listB = self.read_VMI_one_file(root, sufix, which, showmessage=showmessage)
        if not ret:
            return ret, msg, None, None
        
        return ret, msg, listA, listB
       

    def read_VMI_one_file(self, root, sufix, which, showmessage=False):
        filename = root + 'VMI_%s_%s.txt'%(sufix, which)
        
        if not os.path.isfile(filename):
            print 'read_VMI_one_file: error', sufix
            return False, "Could not find: '" + filename +"'", None, None
        try:
            f = open(filename, 'r')
        except IOError, e:
            print 'error', e.errno, e.strerror
            return False, "Could not open: '" + self.completeFilename +"'", None, None
            
        try:            
            vmi_dic = eval(f.read())
            f.close()
            
        # except IOError, e:
        except Exception as e:
            print 'error eval(f.read())', e  # e.errno, e.strerror
            return False, "Could read: '" + self.completeFilename +"'", None, None
    
        if showmessage:
            print '\n-------------------------------' 
            print '  reading keys' 
            keys = vmi_dic.keys()
            
            print 'keys =', keys, ' - ', type(keys)
            print '-------------------------------' 


        if which == 'mi':
            return True, "VMI() crosscor Table Opened " + filename, \
                vmi_dic['MI'], vmi_dic['MICorr']
        elif which == 'se':
            return True, "SE(MI) Table Opened " + filename, \
                vmi_dic['SE'], vmi_dic['SEcorr']
        elif which == 'hShannon':
            return True, "HShannon Table Opened " + filename, \
                vmi_dic['H'], vmi_dic['HCorr']
        elif which == 'seh':
            return True, "SE(h) Table Opened " + filename, \
                vmi_dic['SEH'], vmi_dic['SEHcorr']
        else:
            return True, "Pij Table Opened " + filename, \
                vmi_dic['Pi'], vmi_dic['IJ']
                
                
    def write_HMI_files(self, sufix, MIlist, listMIcorr, SeMIList, SeMIcorrList,listN):
        
        values = dict()
        values['MI'] = tuple(MIlist)
        values['MIcorr'] = tuple(listMIcorr)
        self.write_HMI_one_file(sufix, 'mi', 'txt', values)

        values = dict()
        values['SE'] = tuple(SeMIList)
        values['SEcorr'] = tuple(SeMIcorrList)
        values['N'] = tuple(listN)
        self.write_HMI_one_file(sufix, 'se', 'txt', values)

        '''
        values = dict()
        values['Pij'] = tuple(listPij)
        values['Pi'] = tuple(seqPis)
        self.write_HMI_one_file(sufix, 'pij', 'txt', values)
        '''
        
    def write_HMI_one_file(self, sufix, mode, typeEnd, dic):
        try:
            filename = self.root + 'HMI_%s_%s.%s'%(sufix, mode, typeEnd)
            
            f = open(filename, 'w')
        except:
            print 'Could not save', filename
            return False
        
        print "write_HMI_one_file", filename
        try:
            f.write(str(dic))
            f.flush()  
                                       
            ret = True          
        except:
            print 'Could not save', filename, 'error in str(dic): memory?'
            ret = False
        finally:       
            f.close()
            return ret
                

    '''  sufix, self.dicPiList, self.HShannonList, self.HShannonCorrList,  self.SeHShannonList, self.SeHShannonCorrList, \
         self.MIlist, self.MIcorrList,  self.SeMIList, self.SeMICorrList, self.ijList, showmessage  '''
    def write_VMI_files(self, sufix, dicPiList, HShannonList, HShannonCorrList, SeHShannonList, SeHShannonCorrList, \
                        MIlist, MIcorrList, SeMIList, SeMICorrList, ijList, showmessage=False):
        
        values = dict()
        values['H'] = tuple(HShannonList)
        values['HCorr'] = tuple(HShannonCorrList)
        self.write_VMI_one_file(sufix, 'hShannon', 'txt', values)

        values = dict()
        values['SEH'] = tuple(SeHShannonList)
        values['SEHcorr'] = tuple(SeHShannonCorrList)
        self.write_VMI_one_file(sufix, 'seh', 'txt', values)
        
        values = dict()
        values['MI'] = tuple(MIlist)
        values['MICorr'] = tuple(MIcorrList)
        self.write_VMI_one_file(sufix, 'mi', 'txt', values)

        values = dict()
        values['SE'] = tuple(SeMIList)
        values['SEcorr'] = tuple(SeMICorrList)
        self.write_VMI_one_file(sufix, 'se', 'txt', values)

        values = dict()
        values['Pi'] = tuple(dicPiList)
        values['IJ'] = tuple(ijList)
        self.write_VMI_one_file(sufix, 'pij', 'txt', values)


    def write_VMI_one_file(self, sufix, mode, typeEnd, dic):
        try:
            filename = self.root + 'VMI_%s_%s.%s'%(sufix, mode, typeEnd)
            f = open(filename, 'w')
        except:
            print 'Could not save', filename
            return False
        
        try:
            f.write(str(dic))
            f.flush()  
                                      
            ret = True          
        except:
            print 'Could not save', filename, 'error in str(dic): memory???'
            ret = False
        finally:       
            f.close()
            return ret
        
    def write_file(self, root, filename, stri, showSaved=False):
        filename = root + filename
        
        try:
            f = open(filename, 'w')
            f.write(stri)
            f.flush()
            ret = True
            if showSaved:
                print("File saved: %s"%(filename))
        except:
            print 'Could not write in ' + filename
            ret = False
        finally:
            f.close()
            
        return ret
       
    def read_file(self, root, filename):
        
        path_filename = root + filename

        if not os.path.isfile(path_filename):
            print 'Could not find', path_filename
            return False, None
        try:
            f = open(path_filename, 'r')
        except IOError, e:
            print 'Could not open %s, error no: %i, error: %s'%(path_filename, e.errno, e.strerror)
            return False, None
            
        try:            
            lista = f.read()
            ret = True
        except Exception as e:
            print "Could not read: '%s'"%(path_filename)
            ret = False
            lista = []
        finally:
            f.close()
        
        return ret, lista
    
                       
    def write_file_data_dic(self, data_dict):
        try:
            f = open(self.completeFilename, 'w')
            
            s = str(data_dict)
            f.write(s)
            f.flush()
            ret = True
        except:
            print 'Could not write in ' + self.completeFilename
            ret = False
        finally:
            f.close()
            
        return ret

    def write_file_data_dic2(self, filename, dic):
        try:
            f = open(filename, 'w')
            
            s = str(dic)
            f.write(s)
            f.flush()
            ret = True
        except:
            print 'Could not write in ' + self.completeFilename
            ret = False
        finally:
            f.close()
            
        return ret
    
    def read_file_data_dic(self):
        try:
            f = open(self.completeFilename, 'r')
            result = eval(f.read())
            f.close()
        except:
            return False, "Could not open" + self.completeFilename
    
        return True, result


    def save_tsallis_file(self, numSequences, numOfExperiments, qs, ms, val, sdv):
        values = dict()
        
        dummy =[]
        dummy.append(numSequences)
        values['numSequences'] = tuple(dummy)

        dummy =[]
        dummy.append(numOfExperiments)
        values['numOfExperiments'] = tuple(dummy)
        
        values['q'] = tuple(qs)
        values['m'] = tuple(ms)
            
        for ai in range(len(qs)):
            seqValResult = []
            seqSdvResult = []
        
            for i in range(len(ms)):
                seqValResult.append(val[i][ai]) 
                seqSdvResult.append(sdv[i][ai]) 
                
            values['val'+str(qs[ai])] = seqValResult
            values['sdv'+str(qs[ai])] = seqSdvResult
        
        self.write_file_data_dic(values)

    def read_tsallis_random_file(self, showKeys = False):
        try:
            f = open(self.completeFilename, 'r')
            self.tsallis_dic = eval(f.read())
            f.close()
            
            if showKeys:
                print '\n-------------------------------' 
                print '  reading keys' 
                keys = self.tsallis_dic.keys()
                
                print 'keys =', keys, ' - ', type(keys)
                print '-------------------------------' 
                
                for i in range(len(keys)):
                    key = keys[i]
                    print key 
        except:
            return False, "Could not open Tsallis Random Table " + self.completeFilename
    
        return True, "Ok, Tsallis Random Table opened " + self.completeFilename
    
    def pos_tsallis_random_file(self, posi, q, showmessage = False):
        if not self.tsallis_dic:
            self.read_tsallis_random_file(showKeys=showmessage)
            
        try:
            value = self.tsallis_dic['val'+str(q)][posi-1]
            sdv   =  self.tsallis_dic['sdv'+str(q)][posi-1]
            
            if showmessage:
                print 'posi', posi, 'q=',q, '<h>', value, 'sdv', sdv
        except:
            print 'Position',posi,' or q', q, 'key', 'val'+str(q), 'does not exist'
            return None, None


        return value, sdv

    def save_shannon_file(self, numSequences, numOfExperiments, ms, val, sdv):
        values = dict()
        
        dummy =[]
        dummy.append(numSequences)
        values['numSequences'] = tuple(dummy)

        dummy =[]
        dummy.append(numOfExperiments)
        values['numOfExperiments'] = tuple(dummy)
        
        values['m'] = tuple(ms)
        values['val'] = tuple(val)
        values['sdv'] = tuple(sdv)
        
        self.write_file_data_dic(values)

        
    def pos_shannon_random_file(self, posi, showmessage = False):
        if not self.shannon_dic:
            print 'Shannon file does not exist'
            return None, None
            
        try:
            value = self.shannon_dic['val'][posi-1]
            sdv   = self.shannon_dic['sdv'][posi-1]
            
            if showmessage:
                print 'posi', posi, '<h>=',value
                print 'posi', posi, 'sdv=',sdv
        except:
            print 'Position does not exist'
            # print 'shannon_dic ', self.shannon_dic
            return None, None

        return value, sdv
    
        
    def MI_prefix_sufix_txt(self, prefix_sufix, which):
        filename  = prefix_sufix + '_' + which + '.txt'
        return filename
        
    def MI_name_txt(self, filename, which):
        filename = filename.replase(".fasta","")
        filename = filename.replase(".fas","")
        filename = filename.replase(".txt","")
        
        return filename + '_' + which + '.txt'






''' old class ProteinFunctionNew: 
, numLines, numCols, legColumns, 
                 legendTitle='No legend', title="No Title"'''
class Entropy:
    def __init__(self, desk):

        self.basic = Basic()
        self.numSequences = 0

        self.words = None
        self.doubleWords = None
        
        self.mySequence = MySequence(desk, showmessage=False)
        self.clearEntropySeqs()

        '''
        self.numLines = numLines
        self.numCols = numCols
        self.legColumns = legColumns
        self.legendTitle = legendTitle
        self.title = title    
        '''
        
        if (desk.showmessage):
            print 'Initializing', self.filename
            
    def read_shannon_random_file(self, root, filename, letter, showKeys = False):
        filename = root+filename
        try:
            f = open(filename, 'r')
            
            if letter == 1:
                self.shannon_dic1 = eval(f.read())
            else:
                self.shannon_dic2 = eval(f.read())
                
            f.close()
            
        except:
            print 'error Could not open Shannon Random Table', filename
            return False, "Could not open Shannon Random Table " + filename
    
        if showKeys:
            print '\n-------------------------------' 
            print '  reading keys'
            
            if letter == 1: 
                keys = self.shannon_dic1.keys()
            else:
                keys = self.shannon_dic2.keys()
            
            print 'keys =', keys
            print '-------------------------------' 

        return True, "Ok, Shannon Random Table " + filename + ' opened.'

    def clearEntropySeqs(self):
        self.HShannonList = [] 
        self.entropyTsallis = [] 
        self.entropyRenyi = []  
        self.diversitySimpson = [] 
        
        self.matCounting = [] 
        self.mutualMap = [] 

    # count letter from protein and increment self.matCounting
    def countLetters(self, showmessage):
        showmessage = True
        if (showmessage):
            print 'countLetters ... len ', len(self.mySequence.seqTranspose)
            print 'inner ... len ', len(self.mySequence.seqTranspose[0])
        
        for row in range(len(self.mySequence.seqTranspose)):
            # string has all different letters
            string = ''
            stringList = []
            
            if (showmessage):
                print '--------------------------------------'
            
            for column in range(len(self.mySequence.seqTranspose[row])):
                letter = self.mySequence.seqTranspose[row][column]
                # did y find this letter?
                pos = string.rfind(letter)
                
                if pos < 0:
                    if (showmessage) and (row < 3):
                        print 'ins in ', string, '<-'
                   
                    string = string + letter
                    '''
                      mat4 = row, letter, count, porcentage
                    '''
                    mat4 = [row, letter, 1, 0.00]
                    stringList.append(mat4)
                    
                    if (showmessage):
                        print 'ins row ', mat4[0], " pos ", mat4[1], " letter '", mat4[2], "' count = ", mat4[3]

                else:
                    if (showmessage) and (row < 3):
                        print 'upd in ', string, 'pos', pos, 'len', len(string)

                    mat4 = stringList[pos]
                    
                    mat4[3] += 1
                    if (showmessage):
                        print 'upd row ', mat4[0], " pos ", mat4[1], " letter '", mat4[2], "' count = ", mat4[3]
    
            self.matCounting.append(stringList)


    def countManyLetters(self, seqIndivs, numOfLetters, showmessage):
        lenTrans = len(seqIndivs) / numOfLetters
        
        if (showmessage):
            print '--------------------------------------'
            print 'seqIndivs',seqIndivs
            print 'count polissilables ... len = len(seqIndivs) / numOfLetters =', lenTrans
            print '      inner ... len ', len(seqIndivs[0])
            print '--------------------------------------'
        
        numSequences = len(seqIndivs[0])
        
        for auxRow in range(lenTrans):
            row = auxRow * numOfLetters
            
            # string has all different letters/monosilabs
            string = []
            stringList = []

            for column in range(numSequences):
                letter = ''
                
                # sum letters in vertical
                for pos in range(numOfLetters):
                    letter = letter + seqIndivs[row + pos][column]
                    
                if (showmessage):
                    print 'letter', letter
                    
                # did y find this letter?
                found = False
                for posStri in range(len(string)):
                    if (letter == string[posStri]):
                        found = True
                        break
                
                if not found:
                    string.append(letter)
                    mat4 = [auxRow, letter, 1, 0.00]
                    stringList.append(mat4)

                    if (showmessage):
                        print ' -->*1 ins row =', mat4[0], " letter =", "'" + mat4[1] + "'", " count =", mat4[2], " % =", mat4[3]

                else:
                    mat4 = stringList[posStri]
                    mat4[2] += 1
                    if (showmessage):
                        print ' -->*2 ins row =', mat4[0], " letter =", "'" + mat4[1] + "'", " count =", mat4[2], " % =", mat4[3]

            if (showmessage):
                print ' -->*3 stringlist ', stringList
                
            self.matCounting.append(stringList)
    
    ''' Flavio: 12/02/2014 '''
    def countManyLettersDic(self, seqs, numOfLetters):
        # num_mers = len(seqIndivs) / numOfLetters
        ''' the matrix IS NOT transposed '''
        lcols = len(seqs[0]) - (numOfLetters - 1)
        lseqs = len(seqs)
        
        words = self.build_words(numOfLetters)

        for col in range(lcols):
            dic = {}
            for word in words:
                '''         qtt, % '''
                dic[word] = [0, 0.00]
                                
            for row in range(lseqs):
                word = str(seqs[row])[col:col+numOfLetters]
                ''' search for this word in each line '''
                dic[word][0] += 1

           
            for word in dic.keys():
                dic[word][1] += dic[word][0] / float(lseqs)
                
            self.matCounting.append( [ [key, dic[key][0], dic[key][1]] for key in dic.keys()] )

    def build_words(self, num):
        nuc = ['A','G','T','C']
        
        if num <= 1:
            return nuc
        
        
        seq = copy.deepcopy(nuc)
        
        
        while num <> 1:
            num -= 1
            seq_aux = []
            for nuc1 in seq:
                for nuc2 in nuc:
                    seq_aux.append(nuc1+nuc2)
                    
            seq = seq_aux
            
        return seq
        
    ''' Flavio: 12/02/2014 '''
    def countTwoSiteWords(self, seqs, posI, posJ, numOfLetters):
        # num_mers = len(seqIndivs) / numOfLetters
        ''' the matrix IS NOT transposed '''
        lseqs = len(seqs)
        
        if not self.doubleWords:
            self.doubleWords = self.build_words(numOfLetters*2)
        '''
        for word in self.words:
            print word
        '''

        dic = {}
        for word in self.doubleWords:
            '''          qtt, % '''
            dic[word] = [0, 0.00]
                            
        for row in range(lseqs):
            word = str(seqs[row])[posI:posI+numOfLetters] + str(seqs[row])[posJ:posJ+numOfLetters]
            ''' search for this word in each line '''
            dic[word][0] += 1
            
        ''' could be fast '''

       
        for word in dic.keys():
            dic[word][1] += dic[word][0] / float(lseqs)
            
        return dic



    def calcShannonDoubleVerticalEntropy(self, seqs, posI, posJ, numOfLetters):

        dic = self.countTwoSiteWords(seqs, posI, posJ, numOfLetters) 
        '''
        if posI == 9 and posJ == 371:
            for key in dic.keys():
                print key, dic[key]
                
            print 'yes...'
        '''
        
        listPijs = []
        hShan = 0 
        totVar = 0
        totVar2 = 0

        for key in dic.keys():
            pi = dic[key][1]
            listPijs.append(pi)

            if (pi != 0.) and (pi != 1.):
                hShan -= pi * np.log(pi)
                totVar += pi * (1-pi) * (1. + np.log(pi)) ** 2
                totVar2 += pi * (1-pi)

        return hShan, listPijs, totVar, totVar2


    def calcShannonDoubleVerticalEntropy_WithoutVars(self, seqs, posI, posJ, numOfLetters):

        dic = self.countTwoSiteWords(seqs, posI, posJ, numOfLetters) 

        
        listPijs = []
        hShan = 0 

        for key in dic.keys():
            pi = dic[key][1]
            listPijs.append(pi)

            if (pi != 0.) and (pi != 1.):
                hShan -= pi * np.log(pi)

        return hShan, listPijs
    
        
    def calcVerticalMutualInfo(self, desk, seqs, showmessage=False):
        ''' -------------------------------------------------------
        --- Calculating entropy
        ------------------------------------------------------- '''
        self.desk = desk
        #self.nuc_aa_List = desk.nuc_aa_List
        self.numOfLetters = desk.numOfLetters
        self.isProtein = desk.isProtein
        
        ''' calc entropies e variations '''
        self.calcHSannon_bias_correction(seqs, desk.numOfLetters)
        '''
            self.dicPiList, self.HShannonList, self.HShannonCorrList, \
                            self.SeHShannonList, self.SeHShannonCorrList    
        '''
        MIlist = []
        MIcorrList = []
        ijList = []
        SeMIList = []
        SeMICorrList = []
        
        maxL = len(seqs[0]) - (desk.numOfLetters-1)
        
        desk.showmsg_obs('- Calc join Entropies --- len entropies: %i'%(len(self.HShannonList)))
        
        for i in range(maxL-1):

            if i%50 == 0:
                if desk.showmsg_obs:
                    desk.showmsg_obs(str(i), same_line = True)
                else:
                    print str(i), 

            for j in range(i+1, maxL):

                ''' correction: 09-13/10/2014 '''
                MIk, MIkCorr, seMI, seMICorr = self.crossCorrelation_ij(seqs, i,j, desk.numOfLetters)                

                if MIk < 0:
                    if abs(MIk) <= .000001:
                        MIk = 0.00
                        MIkCorr = 0.00
                        print("*** Correcting *** VMI near zero = %f, is less then 0 in i=%i j=%i"%(MIk, i, j) )
                    else:
                        print("*** Impossible *** VMI = %f, is less then 0 in i=%i j=%i"%(MIk, i, j) )
                        
                    # raw_input('Enter to continue.')
                
                MIlist.append(MIk)
                MIcorrList.append(MIkCorr)
                SeMIList.append(seMI)
                SeMICorrList.append(seMICorr)
                ijList.append([i,j])

        return  self.dicPiList, self.HShannonList, self.HShannonCorrList, \
                self.SeHShannonList, self.SeHShannonCorrList, \
                MIlist, MIcorrList, SeMIList, SeMICorrList, ijList
                
                            
    def fast_calcHSannon(self, X):
        L = float(len(X))

        probs = [X.count(c)/L for c in set(X)]
        
        return np.sum(-p * np.log(p) for p in probs)

    
    def calcHSannon(self, seqs, numOfLetters):
        self.dicPiList = []
        self.HShannonList = []

        lenSeq = len(seqs[0])
        numOfSeqs = len(seqs)
        maxL = lenSeq-numOfLetters+1
        
        ''' for each nucleotide position (residue) '''
        for i in range(maxL):
            hShan = 0
            dicPi = {}
            
            ''' for each row '''
            for row in range(numOfSeqs):
                try:
                    dicPi[seqs[row][i:i+numOfLetters]] += 1
                except:
                    dicPi[seqs[row][i:i+numOfLetters]] = 1
                
            for key in dicPi.keys():
                dicPi[key] /= float(numOfSeqs) 
                ''' varPi += pi * (1-pi) * (np.log(pi) + 1.) ** 2 '''
                
                pi = dicPi[key]
                if (pi != 0.) and (pi != 1.):
                    hShan -= pi * np.log(pi)
                    
                    
            self.dicPiList.append(dicPi)
            self.HShannonList.append(hShan)

            
    def calcHSannon_bias_correction(self, seqs, numOfLetters):
               
        self.dicPiList = []
        self.HShannonList = []
        self.HShannonCorrList = []
        self.SeHShannonList = []
        self.SeHShannonCorrList = []

        lenSeq = len(seqs[0])
        numOfSeqs = len(seqs)
        maxL = lenSeq-numOfLetters+1
        
        for i in range(maxL):
            dicPi = {}
            varPi = 0
            hShan = 0
            
            for row in range(numOfSeqs):
                try:
                    dicPi[seqs[row][i:i+numOfLetters]] += 1
                except:
                    dicPi[seqs[row][i:i+numOfLetters]] = 1
                
            for key in dicPi.keys():
                dicPi[key] /= float(numOfSeqs) 
                
                pi = dicPi[key]
                if (pi != 0.) and (pi != 1.):
                    hShan -= pi * np.log(pi)
                    varPi += (1. + np.log(pi))**2  * pi * (1.-pi)
                    
            B = len(dicPi.keys())
            N = float(numOfSeqs)
            
            ''' Roulston Eq. 13 '''
         
            hShanCorr = hShan + ((B-1)/ (2*N))
            
            VarCorr = 0
            ''' Roulston Eq. 40 '''
            for key in dicPi.keys():
                pi = dicPi[key]
                if (pi != 0.) and (pi != 1.):
                    VarCorr += (np.log(pi)+hShan)**2 * pi * (1-pi)
                
            SeHShannon     = hShan * np.sqrt(varPi/N)
            SeHShannonCorr = np.sqrt(VarCorr/N)
            
            self.dicPiList.append(dicPi)
            self.HShannonList.append(hShan)
            self.HShannonCorrList.append(hShanCorr)
            self.SeHShannonList.append(SeHShannon)
            self.SeHShannonCorrList.append(SeHShannonCorr)
                                    
    def crossCorrelation_ij(self, seqs, i,j, numOfLetters):
        dicPij = {}
        seqsLen = len(seqs)
        N = float(seqsLen)
        
        ''' scan all AGTC or AA 
            Var = 1/n * sum( pij *(1-pij) )

            SUM OF 3 FRAMES '''
        
        for row in range(seqsLen):
            n1 = seqs[row][i:i+numOfLetters] 
            n2 = seqs[row][j:j+numOfLetters] 

            try:                      
                dicPij[n1+n2] += 1
            except:
                dicPij[n1+n2] = 1
            
        for key in dicPij.keys():
            dicPij[key] = dicPij[key] / N  

        '''
            have been calculated at self.calcHSannon(seqs, numOfLetters)
            self.dicPiList.append(dicPi)
        '''
        dicPi = self.dicPiList[i]
        Bi = sum([1 for val in dicPi.values() if val > 0])
             
        dicQj = self.dicPiList[j]
        Bj = sum([1 for val in dicQj.values() if val > 0])
        
        Bij = sum([1 for val in dicPij.values() if val > 0])
        
        
        correction = (Bi+Bj-Bij-1) / (2.*N)
        # print 'Bij %i,  Bi %i  Bj %i, and N = %i, correction=%f'%(Bij, Bi, Bj, N, correction)
        
        totMI = 0.
        totVar  = 0.
        totVarCorr = 0.

        for key in dicPij.keys():
            n1 = key[0:numOfLetters]
            n2 = key[numOfLetters:numOfLetters+numOfLetters]
    
            Pij = dicPij[key]
            Pi = dicPi[n1]
            Qj = dicQj[n2]

            if (Pij > 0 and Pij < 1) :
                    
                mi = Pij * np.log(Pij / (Pi * Qj) )

                ''' Roulston equation 42 se = sqrt(var/N)'''
                totVarCorr += (np.log(Pi) + np.log(Qj) - np.log(Pij) + mi)**2 * (Pij*(1.-Pij))
                totVar     += mi**2 * (Pij*(1.-Pij))
                totMI += mi
    
        ''' 
        if self.HShannonList[i] > 0:
            VarHi = (self.SeHShannonList[i]/self.HShannonList[i])**2
        else:
            VarHi = 0.
            
        if self.HShannonList[j] > 0:
            VarHj = (self.SeHShannonList[j]/self.HShannonList[j])**2
        else:
            VarHj = 0.
            
        totVar = totVarHij + N*( VarHi + VarHj)
        '''
                
        '''   MI cannot be negative never: 09/06/2015  '''
        totMIcorr = totMI + correction
        if totMIcorr < 0:
            # print('Bias correction turn to negative !', totMIcorr, totMI, correction)
            totMIcorr = 0
        
        return totMI, totMIcorr, np.sqrt(totVar/N), np.sqrt(totVarCorr/N)
   
    

    def printAllEntropy(self, iLoop, numOfLetters, randShannonEntropy, randShannonSdv, rootName, wantsTsallisRenyi, output, graph1, graph2):
        self.calcEntropy(seqIndivs=self.mySequence.seqTranspose, numOfLetters=numOfLetters, randShannonEntropy=randShannonEntropy, showmessage=False, wantsShannon=True, wantsTsallisRenyi=wantsTsallisRenyi)
        
        graph1.histogramBar('HIST', 1, 2, 1, title='H Shannon', alfaValue=None, seqAlfa=None, seqName=self.name, seqX=None, seqY=self.HShannonList, numOfLetters=numOfLetters, randAvg=randShannonEntropy, randStd=randShannonSdv, showError=True)
        print 'randEntropy', randShannonSdv
        
        graph1.densityBar(lines=1, columns=2, numOfFig=2, title='H Shannon', seq=self.HShannonList, numOfLetters=numOfLetters, randStd=randShannonSdv, acumuleError=True)
        print 'end shannon'
        

        lines = int(len(self.alfa) / 2.)
        cols  = 2
        
        if wantsTsallisRenyi == 'R':
            for j in range (lines*cols):
                seqAny = []
                for i in range (len(self.HShannonList)):
                    # print j, i, 'self.entropyRenyi[i][j]', self.entropyRenyi[i][j]
                    seqAny.append(self.entropyRenyi[i][j])
                    
                graph2.histogramBar('HIST', lines, cols, j+1, 'H Renyi', self.alfa[j], None, self.name, None, seqAny, numOfLetters, 0, 0, False)
    
        else:    
            for j in range (lines*cols):
                seqAny = []
                
                # print '****', j, 'len alfa', len(self.alfa), 'loop', len(self.HShannonList)
                for i in range (len(self.HShannonList)):
                    # print j, i, 'self.entropyTsallis[i][j]'
                    seqAny.append(self.entropyTsallis[i][j])
                    
                graph2.histogramBar('HIST', lines, cols, j+1, 'H Tsallis', self.alfa[j], None, self.name, None, seqAny, numOfLetters, 0, 0, False)


    
    def printEntropyDistribution(self, iLoop, numLines, numOfLetters, randEntropy, randSdv, showError, root_image, output, graph, oneByOne, CONST):
        self.calcEntropy(self.mySequence.seqTranspose, numOfLetters, randEntropy, False)
    
        entropyCurve = []
        
        for j in range (len(self.alfa)):
            seqAny = []
            for i in range (len(self.HShannonList)):
                seqAny.append(self.entropyRenyi[i][j])
                
            np.mean = CONST.np.mean(seqAny) / randEntropy
            entropyCurve.append(np.mean)
            ''''    
            if (self.alfa[j] == 0.99):
                seqX.append(1)
                np.mean = CONST.np.mean(self.HShannonList)
                entropyCurve.append(np.mean)
            '''       
        
        if oneByOne:
            lines = 1
            columns = 2
            figureNum = 1
        
            graph.histogramBar('HIST', lines, columns, figureNum, 'H Shannon', None, None, self.name, None, self.HShannonList, numOfLetters, randEntropy, randSdv, showError)
            graph.histogramBar('LINE', lines, columns, figureNum+1, 'Renyi Entropy', None, None, self.name, self.alfa, entropyCurve, 1, randEntropy, randSdv, False)
            graph.show_save(self.name, numOfLetters, output,'png', root_image)   
            
        else:
            lines = numLines
            columns = 2
            figureNum = 1+(2*iLoop)
            
            graph.histogramBar('HIST', lines, columns, figureNum, 'H Shannon', None, None, self.name, None, self.HShannonList, numOfLetters, randEntropy, randSdv, showError)
            graph.histogramBar('LINE', lines, columns, figureNum+1, 'Renyi Entropy', None, None, self.name,  self.alfa, entropyCurve, 1, randEntropy, randSdv, False)
        
            if (iLoop == (numLines-1)):
                graph.show_save(self.name, numOfLetters, output,'png', root_image)                     


    # cuidado pode ter melado outra printEntropyTRCurves // tem que mudar de nome
    def printEntropyTRCurves(self, iLoop, numFigure, entTsallisFile, numOfLetters, randShannonEntropy, root, wantsTsallisRenyi, printLegend=True, showmessage=False):
        # for each species iLoop
        self.calcEntropy(seqIndivs=self.mySequence.seqTranspose, numOfLetters=numOfLetters, randShannonEntropy=randShannonEntropy, showmessage=showmessage, wantsShannon=True, wantsTsallisRenyi=wantsTsallisRenyi)
        
        entropyCurve = []
        
        print 'numSequences', self.mySequence.numSequences
        
        for j in range (len(self.alfa)):
            seq = []
            
            q = self.alfa[j]

            # todo !!! Renyi pos_renyi_random_file
            randEntropy, randSdv = entTsallisFile.pos_tsallis_random_file(self.mySequence.numSequences, q, showmessage=False)
            print 'q',q, 'randEntropy', randEntropy, 'randSdv', randSdv
            
            if (not randEntropy) or (not randSdv):
                print 'Couldnt read Tsallis error table !!! for q = ', q
                randEntropy = 0
                randSdv = 0
            
            for i in range (len(self.HShannonList)):

                if (wantsTsallisRenyi == 'R'):
                    seq.append(self.entropyRenyi[i][j])
                else:
                    seq.append(self.entropyTsallis[i][j])
            
            # entropyCurve has np.mean + vc
            if randEntropy > 0:
                value = [self.basic.maximum(seq)*self.basic.convShannonTsallis, randEntropy*self.basic.convShannonTsallis, randSdv*self.basic.convShannonTsallis]
            else:
                value = [self.basic.maximum(seq)*self.basic.convShannonTsallis, 0., 0.]
                
            print ' ** ', value
            entropyCurve.append(value)

        subTitle = ''
        label = self.name.replace("Drosophila ", "")
        label = label.replace("Drosophila_", "")
        label += ' ' + str(self.numSequences) + ' seqs'
        
        xlabel='q'
                 
        if (wantsTsallisRenyi == 'R'):
            ylabel='Hmax Renyi'
        else:
            ylabel='Hmax Tsallis (Sq/k * 1/ln(2))'

        self.graph.multilineGeneral_TR_WithTitle(iLoop=iLoop, numFigure=numFigure, seqX=self.alfa, entropyCurve3Val=entropyCurve, title=subTitle, label=label, printLegend=printLegend, xlabel=xlabel, ylabel=ylabel)

    # cuidado pode ter melado outra printEntropyTRCurves // tem que mudar de nome
    def printEntropyShannonHistogram(self, iLoop, numFigure, winWidth, numOfLetters, randShannonEntropy, randShannonSdv, root, subTitle='', printLegend=True,showmessage=False):
        # for each species iLoop
        
        lenTranspose = len(self.mySequence.seqTranspose)
        lenWidth = len(self.mySequence.seqTranspose[0])
        num_windows= int(lenTranspose / winWidth)
        
        if showmessage:
            print 'num_windows', num_windows, 'lenTranspose', lenTranspose, 'lenWidth',lenWidth
            print self.mySequence.seqTranspose
        
        entropyCurve = []
        seqX=[]

        # getting numOfLetters rows
        for j in range(num_windows):
            seqWindow = []
            # print '---', winWidth
            for i in range(winWidth):
                seqWindow.append(self.mySequence.seqTranspose[ (j*winWidth) + i])
                # print '***', seqWindow[i]
                
            if showmessage:
                print 'seqWindow', seqWindow
            self.calcEntropy(seqIndivs=seqWindow, numOfLetters=numOfLetters, randShannonEntropy=randShannonEntropy, showmessage=showmessage, wantsShannon=True, wantsTsallisRenyi='')
            
            
            entropyCurve.append( [self.basic.np.mean(self.HShannonList), randShannonEntropy, randShannonSdv])
            seqX.append((j+1)* winWidth)
            
            
        subTitle = subTitle
        label = self.name + ' ' + str(lenWidth) + 'indvs'
        
        xlabel=''
        ylabel='<h> Shannon'

        self.graph.multilineGeneralWithTitle(iLoop=iLoop, numFigure=numFigure, seqX=seqX, meanCurve3val=entropyCurve, title=subTitle, label=label, printLegend=printLegend, xlabel=xlabel, ylabel=ylabel)

   
    def printEntropyTRHistogram(self, iLoop, numFigure, entTsallisFile, numOfLetters, randShannonEntropy, root, wantsTsallisRenyi, printLegend=True, showmessage=False):
        # for each species iLoop
        self.calcEntropy(seqIndivs=self.mySequence.seqTranspose, numOfLetters=numOfLetters, randShannonEntropy=randShannonEntropy, showmessage=showmessage, wantsShannon=True, wantsTsallisRenyi=wantsTsallisRenyi)
        
        entropyCurve = []
        
        print 'numSequences', self.mySequence.numSequences
        
        for j in range (len(self.alfa)):
            seq = []
            
            q = self.alfa[j]

            # todo !!! Renyi pos_renyi_random_file
            randEntropy, randSdv = entTsallisFile.pos_tsallis_random_file(self.mySequence.numSequences, q, showmessage=False)
            print 'q',q, 'randEntropy', randEntropy, 'randSdv', randSdv
            
            if (not randEntropy) or (not randSdv):
                print 'Couldnt read Tsallis error table !!! for q = ', q
                randEntropy = 0
                randSdv = 0
            
            for i in range (len(self.HShannonList)):

                if (wantsTsallisRenyi == 'R'):
                    seq.append(self.entropyRenyi[i][j])
                else:
                    seq.append(self.entropyTsallis[i][j])
            
            # entropyCurve has np.mean + vc
            if randEntropy > 0:
                value = [self.basic.maximum(seq), randEntropy, randSdv]
            else:
                value = [self.basic.maximum(seq), 0., 0.]
                
            print ' ** ', value
            entropyCurve.append(value)

        subTitle = ''
        label = self.name.replace("Drosophila ", "")
        label = label.replace("Drosophila_", "")
        label += ' ' + str(self.numSequences) + ' seqs'
        
        xlabel='q'
                 
        if (wantsTsallisRenyi == 'R'):
            ylabel='Hmax Renyi'
        else:
            ylabel='Hmax Tsallis'

        self.graph.multilineGeneralWithTitle(iLoop=iLoop, numFigure=numFigure, seqX=self.alfa, meanCurve3val=entropyCurve, title=subTitle, label=label, printLegend=printLegend, xlabel=xlabel, ylabel=ylabel)



    def printEntropyTRHistogram_Q(self, iLoop, numFigure, seqX, q, entTsallisFile, numOfLetters, randShannonEntropy, randShannonSdv, root, wantsTsallisRenyi, printLegend=True, showmessage=False):
        # for each species iLoop
        self.calcEntropy(seqIndivs=self.mySequence.seqTranspose, numOfLetters=numOfLetters, randShannonEntropy=randShannonEntropy, showmessage=showmessage, wantsShannon=True, wantsTsallisRenyi=wantsTsallisRenyi)
        
        entropyCurve = []
        
        print 'numSequences', self.mySequence.numSequences
        
        for j in range (len(self.alfa)):
            seq = []
            
            qX = self.alfa[j]

            if qX == q:
                # todo !!! Renyi pos_renyi_random_file
                randEntropy, randSdv = entTsallisFile.pos_tsallis_random_file(self.mySequence.numSequences, q, showmessage=False)
                print 'q', q, 'randEntropy', randEntropy, 'randSdv', randSdv, 'ShannonSdv', randShannonSdv
                
                if (not randEntropy) or (not randSdv):
                    print 'Couldnt read Tsallis error table !!! for q = ', q
                    randEntropy = 0
                    randSdv = 0
                
                for i in range (len(self.HShannonList)):
    
                    if (wantsTsallisRenyi == 'R'):
                        seq.append(self.entropyRenyi[i][j])
                    else:
                        seq.append(self.entropyTsallis[i][j])
                
                    if randEntropy > 0:
                        value = [seq[i]*self.basic.convShannonTsallis, randEntropy*self.basic.convShannonTsallis, randSdv*self.basic.convShannonTsallis]
                    else:
                        value = [seq[i]*self.basic.convShannonTsallis, 0., 0.]
                        
                    entropyCurve.append(value)

        subTitle = "q = " + str(q)
        label = self.name + ' ' + str(self.numSequences)
        
        if numFigure > 4:
            xlabel='m sequences'
        else:
            xlabel=''
                 
        if (wantsTsallisRenyi == 'R'):
            ylabel='Hpos Renyi'
        else:
            ylabel='Sq/k * 1/ln(2) Tsallis'

        self.graph.multilineGeneralWithTitle(iLoop=iLoop, numFigure=numFigure, seqX=seqX, meanCurve3val=entropyCurve, title=subTitle, label=label, printLegend=printLegend, xlabel=xlabel, ylabel=ylabel)


            
    def calc_anova(self, seqMI, sufix):
            
        #one way anova
        f_value, p_value = f_oneway(*seqMI)
        #f_value, p_value = 0,0
        
        stri = sufix + '\n'
        
        if p_value <= 0.05:
            stri += 'At least one distribution is statistically different from the others, by ANOVA'
        else:
            stri +=  'The distributions are statistically similar, by ANOVA'
            
        stri += '\nf_value %f   p_value %.5e' %(f_value, p_value)
        
        return stri + '\n\n'
    
    def print_data_summary(self, desk, speciesParams, roundVal=4, filename=None, stri = ''):
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
        if desk.saveData:
            self.write_summary_file(desk, filename, str_param)

    def write_summary_file(self, desk, filename, stri):
        filename = desk.rootTable + filename
        try:
            f = open(filename, 'w')
            s = str(stri)
            f.write(s)
            f.flush()
            
            desk.showmsg_obs('write summary file %s'%filename)
            
        except:
            desk.showmsg_obs('Could not write param file %s'%filename)
            return False
        
        finally:
            f.close()
        
        
        return True
