'''
Created on 02/09/2013
Updated on 09/07/2014

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica
'''
from classes import GeneBankParser as GBP

#--- problems with PyInstaller
# from PIL import Image
# _ = Image.Image()

class Pipe():
    def __init__(self, desk):

        self.failure = True
        self.dic_species = {}
        self.error_msg = ''
        self.desk = desk
        
        try:
            desk.get_params()
            desk.showGraph = False
        except:
            self.error_msg = 'Could not get parameters.'
            return
        

        if desk.organism=='':
            desk.showmsg_obs('Define at least one Organism.')
            return

        if desk.title=='' and desk.gene=='':
            desk.showmsg_obs('Define at least a gene or a title.')
            return
        
        
        gbp = GBP.GeneBankParser(desk)
        
        if not gbp.read_and_parse(desk=desk, cntShow=5, maxiSeq=desk.maxiSeq, minSeq=desk.minSeq, stop=desk.stop, showmessage=desk.showmessage):
            self.error_msg = 'Error while parsing.'
            return
        
        gbp.show_data(cntShow=5, stop=desk.stop, showmessage=desk.showmessage)
    
        '''
        if len(myExon_Exon) >= cutoffLength and (is_partial == False or includePartial):
        
        minimu length to accept a sequence
        
        in NCBI many stored sequences are partial, others 'total'
        you chose: includePartial
        '''
    
        desk.dicParams = gbp.save_data_cds_exon_intron_prot(desk, includePartial=True, showmessage=False)    
        self.failure = False



            