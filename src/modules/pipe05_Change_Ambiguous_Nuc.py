'''
Created on 13/06/2013
Updated on 02/06/2014

@author: Flavio Lichtenstein

'''
import classes.Select_Species as SSpec
''' regular expression '''
import re
#from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC

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
        ''' for euchariots '''
        startCodon = 'ATG'
        stopCodonList = ['TAG','TAA', 'TGA']
        filterStartStopCodon = desk.start_stop_codon

        listConsensus = []
        iLoop = 0
        
        for minmax in ["mincut", "maxmer"]:
            iLoop += 1

            sSpec = SSpec.Select_Species(desk)

            desk.output_filename_purged =  ('%s_%s_%s_%s_%iL_cutoff%i_purged.fasta') % \
                         (desk.organism, minmax, desk.seqType, desk.gene_title, desk.cutoffLength, desk.cutoffNumSeq)

            filename = desk.output_filename_purged
            
            stri = '%i/2) %s'%(iLoop, filename)
            desk.showmsg_obs(stri)
            
            countGood, countNoSart, countNoStop, listPolymorphism = 0, 0, 0, []
            
            if not sSpec.read(filename, showmessage=desk.showmessage):
                stri = 'Could not find %s. Impossible to go further.'%(filename)
                stri += "\n\nProbably you didn't purged (and aligned?) this file. Run purging method."
                self.error_msg = stri
                print self.error_msg
                return

            
            j = 0

            for seq in sSpec.mySeq.seqs:
                # desc = sSpec.mySeq.idList[j]
                j += 1
                
                if filterStartStopCodon:
                    ''' supposed any frame '''
                    start = [m for m in re.finditer(startCodon, seq)]
                    if start:
                        start = min(start)
                    else:
                        start = -1
                    
                    ''' if in an intron is not an StopCodon ... review '''
                    stop = []
                    for stopCodon in stopCodonList:
                        stop += [m.start() for m in re.finditer(stopCodon, seq)]
                        
                    if stop:
                        stop = min(stop)
                    else:
                        stop = -1
                    
                    good = True
                    if start < 0:
                        if desk.showmessage:
                            print ', there is no start codon',
                        countNoSart += 1
                        good = False
                        
                    if stop < 0:
                        if desk.showmessage:
                            print ',  there is no stop codon',
                        countNoStop += 1
                        good = False
        
                    if good:
                        if desk.showmessage:
                            print '- L=%i, start=%i stop=%s'%(len(seq), start, str(stop)),
                        countGood += 1
                        
                
                whichIUPACpos = [posi for posi in range(len(seq)) if seq[posi] not in ['A','G','T','C']]
                
                if whichIUPACpos:
                    # desk.showmsg_obs('  >>> %s - IUPAC codes'%(desc))
                    for pos in whichIUPACpos:
                        # print pos+1,
                        if pos not in listPolymorphism:
                            listPolymorphism.append(pos)
                    # print ''
                #else:
                    # desk.showmsg_obs('  ok: %s - IUPAC codes'%(desc))
                    
        
            #print ''
            if filterStartStopCodon:
                print '  ## %i good seqs,  %i no start, %i no stop'%(countGood, countNoSart, countNoStop)
        
            if listPolymorphism:
                listPolymorphism.sort()
                # print '     >> Finding Consensus: mismatch list:', listPolymorphism
                listConsensus = sSpec.mySeq.build_DNA_consensus(sSpec.mySeq.seqs, listPolymorphism, force=True)
                
                #hasConsensus = True
                for seq in listConsensus:
      
                    whichIUPACpos = [posi for posi in range(len(seq)) if seq[posi] not in ['A','G','T','C']]
                    if whichIUPACpos:
                        #hasConsensus = False
                        desk.showmsg_obs('     !Impossible! Problems in consensus:', same_line=True)
                        for pos in whichIUPACpos:
                            desk.showmsg_obs(str(pos), same_line=True)
                        desk.showmsg_obs('')

                '''
                if hasConsensus:
                    desk.showmsg_obs('>>> Consensus achieved.')
                '''

            else:
                print '\n>>> No mismatch.\n'
                listConsensus = sSpec.mySeq.seqs
        
            print '-------------- end %s -----------------'%(minmax)
            
               
                            
            seq_records = []
            j = 0
            for seq in listConsensus:
                seq2 = SeqRecord(seq, id=sSpec.mySeq.idList[j]) #, name=name, description=description)
                seq_records.append(seq2)
                j += 1
                
            filename = desk.rootFasta + '%s_%s_%s_%s_%iL_cutoff%i_consensus.fasta'%\
                 (desk.organism, minmax, desk.seqType, desk.gene_title, desk.cutoffLength, desk.cutoffNumSeq)
                        
            sSpec.mySeq.save_fasta_new(filename, seq_records)
    

            if minmax == 'maxmer':
                dic_seqs_maxmer = desk.recalc(seq_records)
            else:
                dic_seqs_mincut = desk.recalc(seq_records)
                
                
            dicSpecies, lista = sSpec.mySeq.split_species_seqs(desk.organism, desk.gene_title, seq_records, showmessage=False)
    
            for species in lista:
                filename = desk.rootFasta + '%s_%s_%s_%s_%s_%iL_cutoff%i_consensus.fasta'%\
                 (desk.organism, minmax, species, desk.seqType, desk.gene_title, desk.cutoffLength, desk.cutoffNumSeq)
                 
                seq_records = []
                for seq_rec, seqId in dicSpecies[species]:
                    seq2 = SeqRecord(seq_rec, id=seqId) #, name=name, description=description)
                    seq_records.append(seq2)
                    
                sSpec.mySeq.save_fasta_new(filename, seq_records)

    
        desk.redefine_params(dic_seqs_mincut, dic_seqs_maxmer, 'consensus')
        
        self.failure = False

