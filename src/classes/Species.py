'''
Created on 06/03/2012
Revised on 20/06/2012

@author: Flavio
'''

import Crud as cCrud
   
class Species:
    def __init__(self):
        
        self.species = []
        self.acc = None
        
    def connect_MySql(self):
        try:
            self.acc = cCrud.Crud()
            self.acc.crud_connect()
            
        except:
            self.acc.crud_close()
            return False

        return True

    def getSpecies(self, mySeq, names):  # featGene
        if (not self.acc):
            return self.species
                
        for i in range(len(names)):
            geneFeature = None
            self.species.append(self.readWallByName(self.acc, mySeq, names[i], geneFeature))
    
        return self.species
    
    
    def copySequence(self, seqOri, seqDestine, iniRow=0, endRow=0):
        if iniRow==0:
            for i in range(len(seqOri)):
                seqDestine.append(seqOri[i])
        else:
            for i in range(endRow - iniRow):
                seqDestine.append(seqOri[iniRow - 1 + i])
                
        return seqDestine
            
    def sumSequence(self, seqOri, seqToSum):
        for i in range(len(seqToSum)):
            seqOri.append(seqToSum[i])
                
        return seqOri
                
    def readWallByName(self, acc, mySeq, name, featGene):
        
        query = "SELECT indv_id FROM Species_Individuals where upper(indv_source) = '" + name.upper() + "'"
        
        '''
        select spec.name, spei.species, indv_id, indv_description, indv_organism organism, 
               gene_searched, feature_gene gene, feature_product product, feature_protein_id protein_id,
               feature_translation seq_protein, sequence seq_dna
        from species_individuals spei
        join species spec on spec.id = spei.species
        where gene_searched like 'ADH%' and
              feature_gene <> 'adh-1' and
              feature_gene <> 'adh1' and
              feature_gene <> 'adh-2' and
              feature_gene <> 'adh2' and
              feature_gene <> 'Adh-Finnegan' and
              feature_gene <> 'adhr' and
              feature_gene <> 'adh-r'
              
        order by spec.name;
        
        '''
        if (featGene):
            query += query + " and feature_gene = '" + featGene + "'";
         
        ret, specieRows = acc.crud_query(query)
        
        if not ret:
            print 'Impossible to read Sql:', query
            return []
        
        return specieRows
            
       
        

        

        
        