'''
Created on 07/08/2013

@author: Flavio Lichtenstein
@local: Unifesp DIS - BioinformaticaSpecie

'''
import os
from BioPythonClass import MySequence
from Bio.SeqRecord import SeqRecord

class Select_Species:
    def __init__(self, desk):
        
            
        self.organism = desk.organism_var.get()
        self.species = desk.species_var.get()
        self.seqType = desk.seqType_var.get()
        self.gene = desk.gene_var.get()
        
        self.cutoffNumSeq = desk.cutoffNumSeq_var.get()
        self.cutoffLength = desk.cutoffLength_var.get()

        self.name = '%s_%s' % (self.organism, self.gene)
        
        self.isProtein = desk.isProtein
        self.dna_prot = desk.dna_prot_var.get()

        self.confirm_gene = desk.confirm_gene_var.get()

        self.rootFasta = desk.rootFasta_var.get()
        self.rootImage = desk.rootImage_var.get()
        self.rootTable = desk.rootTable_var.get()
        
        self.mySeq = MySequence(desk, showmessage=False)


    def read(self, filename, showmessage=False):
        names = []
        names.append([filename, self.name])
        self.maxFig = len(names)

        if not self.mySeq.read_fasta(self.rootFasta + filename):
            return False

        return True
    

    def strange_nucleotide_DNA(self, seq, showmessage=False):
        seq = str(seq)
        lista = []
        for i in range(len(seq)):
            if seq[i] not in ['A','G','T','C']: 
                lista.append([seq[i], i])
            
        if lista:
            return True, lista
        
        return False, lista
    
    def purge_DNA_repeated(self, desk, allocationSymbols=True):
        seq_record_CDS = []
        dicSpecies = {}
        
        listBadTerms=desk.listBadTerms.split(",") 
        showmessage=desk.showmessage
        
        for numSeq in range(self.mySeq.getNumSequences()):
            seq = self.mySeq.seqs[numSeq]
            desc = str(self.mySeq.idList[numSeq])
            
            flag, lista = self.strange_nucleotide_DNA(seq)
            
            ''' once a bad char but in IUPAC - entropy must be calculated by the species consensus '''
            if flag and allocationSymbols:
                flag = False
                for nuc, pos in lista:
                    if nuc not in ['R','Y','W','S','M','K','H','B','V','D','N','-']:
                        flag = True
                        break
                    
            if flag:
                desk.showmsg_obs('!!! %i) Mismatch: %s'%(numSeq, desc))
                for letter, pos in lista:
                    print "   '%s' at %i"%(letter, pos)
                print '   --------------------------'
            else:
                desk.showmsg_obs('>>> %i) No Mismatch: %s'%(numSeq, desc))
                
            
            if showmessage:
                print '   -->', desc
        
            pos = desc.find('_')
            off = 1
            if pos < 0:
                print 'bad structure - description', desc
                return False, None
                              
            seqId = desc[:pos]
            stri = desc[pos+off:]

            pos = stri.find('_')
            if pos < 0:
                print 'bad structure - accession number'
                return False, None

            off = 1
            ''' see http://www.ncbi.nlm.nih.gov/refseq/about/ '''
            if pos < 5:
                ncbi_prefix = stri[pos+off:]
                pos += ncbi_prefix.find('_') + 1
                
            accNum = stri[:pos]
            stri = stri[pos+off:]
            # print seqId, accNum

            repeat = False
            for seqAux in seq_record_CDS:
                if id == seqAux.id:
                    print '  >>> Id repeated', id
                    repeat = True
                    break
                    
                if accNum in seqAux.description:
                    print '  >>> accNum repeated', accNum
                    repeat = True
                    break
                
            ''' repeated are not saved '''
            if repeat:
                continue

            for badTerm in listBadTerms:
                if stri.find(badTerm) > 0:
                    print '  >>> in %s bad term %s found in %s'%(seqId, badTerm, stri)
                    repeat = True
                    break

            ''' repeated are not saved '''
            if repeat:
                continue

            pos = stri.find('|')
            off = 1
            if pos < 0:
                print '  >>> organism: in %s bad structure %s'%(seqId, stri)
                continue
                
            organism = stri[:pos]
            stri = stri[pos+off:]

            pos = stri.find('|')
            off = 1
            if pos < 0:
                print '  >>> species: in %s bad structure %s'%(seqId, stri)
                continue
                
            species = stri[:pos]
            if species == 'arizonae':
                pass
            stri = stri[pos+off:]
            
            off = 1
            pos = stri.find('|')
            if pos < 0:
                print '  >>> gene: in %s bad structure %s'%(seqId, stri)
                continue
                
            gene2 = stri[:pos]
            stri = stri[pos+off:]

            if self.confirm_gene:
                if gene2.upper().strip() <> self.gene.upper().strip():
                    print "in %s %i, gene not confirmed: %s versus %s "%(seqId, numSeq, gene2, self.gene)
                    continue
            
            pos = desc.find('Prod|')
            off = 5
            if pos < 0:
                product = stri
            else:
                product = stri[off:]
            
            if showmessage:
                print species, numSeq, desc
        
            seqId = seqId + \
               '_' + accNum + '_' + organism + '|' + species + '|' + gene2 + '|Prod|' + product
            
            seqrec = SeqRecord(seq, id=seqId) #, name=name, description=description)
            seq_record_CDS.append(seqrec)

            if species in dicSpecies.keys():            
                dicSpecies[species].append(seqrec)
            else:
                dicSpecies[species] = [seqrec]

        return True, dicSpecies

    # def cutoff_seqs(self, output_filename):
    def save_purged(self, output_filename, dicSpecies):
        seqs = []
        
        ''' dicSpecies was prepared at purge_DNA_repeated'''
        for key in sorted(dicSpecies.keys()):
            for seqrec in dicSpecies[key]:
                seqs.append(seqrec)

        self.mySeq.save_fasta_new(output_filename, seqs)
        
        return seqs
    
    def write_file(self, filename, text):
        # saving the file           
        try:
            handle = os.open(filename, "w")
            handle.write(text)
            print "File '" + filename + "' saved."
            ret = True
            
        except ValueError:
            print "File '" + filename + "' not saved. Writing error: ", ValueError
            ret = False

        finally:    
            handle.close()

        return ret

    def read_file(self, filename):
        if not os.path.exists(filename):
            print 'File', filename, 'does not exists.'
            return None
        
        try:
            handle = os.open(filename, "r")
            lines = handle.read()
            print "Reading file '" + filename + "'"
            
        except ValueError:
            print "File '" + filename + "' not found error: ", ValueError
            return None

        finally:    
            handle.close()

        return lines
    
