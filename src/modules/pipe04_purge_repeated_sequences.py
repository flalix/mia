'''
Created on 30/08/2013
Updated on 03/04/2014

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica

'''
from classes import Select_Species as SSpec

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

            
        iLoop = 0
        for minmax in ["mincut", "maxmer"]:
            iLoop += 1

            sSpec = SSpec.Select_Species(desk)

            filenameAlign = ('%s_%s_%s_%s_%iL_cutoff%i_aligned.fasta') %\
                 (desk.organism, minmax, desk.seqType, desk.gene_title, desk.cutoffLength, desk.cutoffNumSeq)
                 
            stri = '%i/2) %s'%(iLoop, filenameAlign)
            desk.showmsg_obs(stri)      
                        
            if not sSpec.read(filename = filenameAlign, showmessage=desk.showmessage):
                stri = '\nCould not find %s. Impossible to go further.'%(filenameAlign)
                stri += "\n\nProbably you didn't aligned this file and renamed as '_aligned'"
                self.error_msg = stri
                print self.error_msg
                return

            ''' creates sSpec.dicSpecies '''    
            ret, dicSpecies = sSpec.purge_DNA_repeated(desk=desk, allocationSymbols=True)
            
            if not ret:
                self.error_msg = 'Error in purge_DNA_repeated() procedure'
                return

            ''' take dic and transf in seqs '''
            filename = desk.rootFasta + filenameAlign.replace("_aligned","_purged")
            seqs = sSpec.save_purged(filename, dicSpecies)
            
            if not seqs:
                self.error_msg = 'No sequences were selected for %s'%(desk.output_filename_purged)
                return

            if minmax == 'maxmer':
                dic_seqs_maxmer = desk.recalc(seqs)
            else:
                dic_seqs_mincut = desk.recalc(seqs)
            
    
        desk.redefine_params(dic_seqs_mincut, dic_seqs_maxmer, 'purge')
        self.failure = False