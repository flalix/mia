#-*- coding: utf-8 -*-
'''
Created on 18/09/2015
Created on 25/09/2015

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica
'''
# import FileDialog
import classes.Classifier as Classifier

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
            
        desk.which = "HMI"
        self.which = "HMI"
            

        desk.unit = 'mnat' if desk.mnat else 'nat'
        desk.factor = 1000 if desk.mnat else 1
        desk.roundVal = 4

        cl = Classifier.Classifier(desk)

        if desk.organism == '':
            self.error_msg = "Define the organism or write 'any'"
            return
        
        if desk.gene == '' and desk.title == '':
            self.error_msg = 'Define at least one Gene or Title'
            return

        ''' -----------------------------------------'''            
        
        ''' can be 0 any k, 1 for frame 1, a
            nd > 1 displaying all frames and 0 '''


        if not desk.de_novo:
            tudo_ok, listH0_JSD, listHa_JSD = cl.read_2lists()
            
            if tudo_ok:
                cl.build_distributions_ROC(listH0_JSD, listHa_JSD)
                
                self.failure = False
                self.error_msg = ''                
                return  
     
        dics = {}
        dics["mincut"] = []
        dics["maxmer"] = []


        #----- fix here ------------
        desk.simulation = False
        #---------------------------

        cl.all_species = {}
        
        for species in desk.dicParams.keys():
            mat = desk.dicParams[species]

            if cl.minmax == "mincut":
                N = mat[2]
            else:
                N = mat[3]
            
            '''
                [N, [], True]
                num of seqs, seqSels, known species
            '''
            if N >= desk.cutoffNumSeq:
                cl.all_species[species] = [N, [], True]

            elif N >= desk.NminSamples:
                cl.all_species[species] = [N, [], False]
               
        
        cl.Lmincut = -1
        cl.Lmaxmer = -1 
        
        if cl.minmax == "mincut":
            cl.LH0 = mat[5]
        else:
            cl.LH0 = mat[6]
            
        cl.Lmincut = mat[5]
        cl.Lmaxmer = mat[6]

        desk.showmsg_obs("Lmincut: %i Lmaxmer: %i bp"%(cl.Lmincut, cl.Lmaxmer) )   

        desk.showmsg_obs(">>>>", cl.minmax )
        listH0_JSD, listHa_JSD = cl.hyp_species_loop(nLoops=desk.numSimLoops)
    
        desk.showmsg_obs("" )
    
        if desk.saveData:
            cl.save_2lists(listH0_JSD, listHa_JSD)
        
        cl.build_distributions_ROC(listH0_JSD, listHa_JSD)

        
        self.failure = False
        self.error_msg = ''
        return



    