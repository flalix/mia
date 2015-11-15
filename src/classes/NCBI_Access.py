'''
Created on 19/08/2012
@author: Flavio Lichtenstein
@local: UNIFESP - bioinformatica

like drosophilaSaveFromNCBI.py
'''

from Bio import Entrez, SeqIO
import BioPythonClass as bioClass
import Crud as cCrud


class NCBI_Access:
    def __init__(self, db, idList):
        self.idList = idList 
        self.db     = db    
        
        self.basic = bioClass.Basic()
        
        Entrez.email = "flalix@gmail.com"
        print 'NCBI_Acess instantiated. db=', idList

    def print_ID_list(self, printResume=True, printDoc=False, printDetails=False, showComment=False):
        print 'printDetails', printDetails
        for i in range( len(self.idList) ):
            
            self.access_ID = self.idList[i]
            
            print "--------------------------------------------------------"
            print 'self.access_ID', self.access_ID
            
            try:
                handle = Entrez.efetch(db=self.db, rettype="gb", id=self.access_ID, retmode="text")
                # print 'handle ', str(handle)
                
                
            except:
                print '   !!!!! error in self.access_ID   !!!!! Stop parsing !!!!', self.access_ID
                continue
        
            if printDoc:
                print "--------------------------------------------------------"
                print handle.read()
                print "--------------------------------------------------------\n"
                handle.close
                # reading again ..
                handle = Entrez.efetch(db=self.db, rettype="gb", id=self.access_ID, retmode="text")

            offset = '  '
            
            try:
                # print 'Parsing NCBI seq ...'
                # print SeqIO.parse(handle, "gb")
                
                # for seq_record in Entrez.parse(handle):
                for seq_record in SeqIO.parse(handle, "gb"):
                    # print 'seq_record....'
                    print "### Parsing '%s', length %i, with %i features" % (seq_record.name, len(seq_record), len(seq_record.features))
        
                    indv_description = seq_record.description
                    indv_name = seq_record.name
                    indv_id = seq_record.id
                    seq_dna = seq_record.seq

                    version = 0
                    str_source = ''
                    comment = ''

                    list_taxonomy = []
                    list_keywords = []
                    list_accessions = []

                    str_date_file_division = ''
                    str_date = ''
                    str_organism = ''
                    str_gi = ''

                    feature_strain = ''
                    feature_mol_type = ''
                    feature_organism = ''
                    feature_db_xref = ''
                    list_genes = ''
                    list_products = ''
                    feature_codon_start = ''
                    seq_protein = ''   # is a list !!
                    list_protein_ids = ''
                    feature_lab_host = ''
                    
                    feature_clone = ''
                    feature_note = ''
                    feature_db_xref = ''
                    feature_tissue_type = ''
                    feature_clone_lib = ''
                            
                    list1_location = []
                    ref1_authors = None
                    ref1_consortium = None
                    ref1_journal = None
                    ref1_medline = None
                    ref1_pubmed = None
                    ref1_comment = None
                    
                    list2_location = []
                    ref2_authors = None
                    ref2_consortium  = None
                    ref2_journal = None
                    ref2_medline = None
                    ref2_pubmed = None
                    ref2_comment = None   
                               
                              
                    i = 1
                    if (printDetails):
                        print 'Keys', seq_record.annotations.keys()
                    for annot in seq_record.annotations:
            
                        j = 1
                        if (annot == 'references'):
                            if (printDetails):
                                print '\n=== Refereces ========================'
                                            
                            annotNumber = 1
                            refObjects = seq_record.annotations[annot]
            
                            for ref in refObjects:
                                # print ref
                                
                                if (annotNumber == 1):
                                        list1_location = ref.location
                                        ref1_authors = ref.authors
                                        ref1_consortium = ref.consrtm
                                        ref1_journal = ref.journal
                                        ref1_medline = ref.medline_id
                                        ref1_pubmed = ref.pubmed_id
                                        ref1_comment = ref.comment
                                        
                                        if (printDetails):
                                            print (offset*(i+1)), str(j) + '. reference1:'
                                            if list1_location:
                                                print (offset*(i+2)), 'location:', list1_location
                                            if ref1_authors:
                                                print (offset*(i+2)), 'authors:', ref1_authors
                                            if ref1_consortium:
                                                print (offset*(i+2)), 'consortium:', ref1_consortium
                                            if ref1_journal:
                                                print (offset*(i+2)), 'journal:', ref1_journal
                                            if ref1_medline:
                                                print (offset*(i+2)), 'medline:', ref1_medline
                                            if ref1_pubmed:
                                                print (offset*(i+2)), 'pubmed:', ref1_pubmed
                                            if ref1_comment:
                                                print (offset*(i+2)), 'comment:', ref1_comment 
                                            
                                            print '------- end ref1 ---------------\n'
                                                                    
                                elif (annotNumber == 2):        
                                        list2_location = ref.location
                                        ref2_authors = ref.authors
                                        ref2_consortium = ref.consrtm
                                        ref2_journal = ref.journal
                                        ref2_medline = ref.medline_id
                                        ref2_pubmed = ref.pubmed_id
                                        ref2_comment = ref.comment
                                        
                                        if (printDetails):
                                            print (offset*(i+1)), str(j) + '. reference2:'
                                            if list2_location:
                                                print (offset*(i+2)), 'location:', list2_location
                                            if ref2_authors:
                                                print (offset*(i+2)), 'authors:', ref2_authors
                                            if ref2_consortium:
                                                print (offset*(i+2)), 'consortium:', ref2_consortium
                                            if ref2_journal:
                                                print (offset*(i+2)), 'journal:', ref2_journal
                                            if ref2_medline:
                                                print (offset*(i+2)), 'medline:', ref2_medline
                                            if ref2_pubmed:
                                                print (offset*(i+2)), 'pubmed:', ref2_pubmed
                                            if ref2_comment:
                                                print (offset*(i+2)), 'comment:', ref2_comment 

                                            print '------- end ref2 ---------------\n'
                                        
                                annotNumber += 1
                                j += 1
                        else:
                            if (printDetails):
                                print '\n==== Annotations ======================='
                            
                            varAny = seq_record.annotations[annot]
                            
                            if (annot == 'comment'):
                                comment = varAny
                                if showComment:
                                    print 'comment ', version
                                
                            elif (annot == 'sequence_version'):
                                version = varAny
                                if (printDetails):
                                    print 'sequence_version ', version, type(version)
                            elif (annot == 'source'):
                                str_source = varAny
                                if (printDetails):
                                    print offset, 'source', varAny
                            elif (annot == 'taxonomy'):
                                list_taxonomy = varAny
                                if (printDetails):
                                    print offset, 'taxonomy', varAny
                            elif (annot == 'date'):
                                str_date = varAny
                                if (printDetails):
                                    print offset, 'date', varAny
                            elif (annot == 'organism'):
                                str_organism = varAny
                                if (printDetails):
                                    print offset, 'organism', varAny
                            elif (annot == 'gi'):
                                str_gi = varAny
                                if (printDetails):
                                    print offset, 'gi', varAny
                            elif (annot == 'accessions'):
                                list_accessions = varAny
                                if (printDetails):
                                    print offset, 'accessions', varAny
                            elif (annot == 'data_file_division'):
                                str_date_file_division = varAny
                                if (printDetails):
                                    print offset, 'data_file_division', varAny
                            else:
                                if (printDetails):
                                    print offset, annot, varAny
                                
                    if (printDetails):    
                        print '\n----------------------------'
                        print "%i features" % len(seq_record.features)
                        # print "%i sub-features," % len(seq_record.SubFeatures)
                      
            
                    for feature in seq_record.features:
                        if (printDetails): 
                            print '\n', offset, str(i) + '. feature:'
                        i += 1
                        
                        j = 2
                        count = 1
                        if (feature.type):
                            stri = feature.qualifiers.get(feature.type)
                            if (stri != None):
                                if (printDetails): 
                                    print (offset*j), count,') type', feature.type, stri
                                count += 1
            
                        if (feature.location):
                            stri = feature.qualifiers.get(feature.location)
                            if (stri != None):
                                if (printDetails): 
                                    print (offset*j), count,') location',  stri
                                count += 1
            
                        if (feature.strand):
                            stri = feature.qualifiers.get(feature.strand)
                            if (stri != None):
                                if (printDetails): 
                                    print (offset*j), count,') strand', feature.strand, stri
                                count += 1
            
                        if (feature.qualifiers):
                            if (printDetails): 
                                print (offset*j), count,') qualifiers:'

                            count += 1
            
                            k = 3
                            count2 = 0
                            for quali in feature.qualifiers:
                                stri = feature.qualifiers.get(quali)
                                count2 += 1
                               
                                if (quali == 'strain'):
                                    feature_strain = stri
                                elif (quali == 'mol_type'):
                                    feature_mol_type = stri
                                elif (quali == 'organism'):
                                    feature_organism = stri
                                elif (quali == 'db_xref'):
                                    feature_db_xref = stri
                                elif (quali == 'gene'):
                                    list_genes = stri
                                elif (quali == 'product'):
                                    list_products = stri
                                elif (quali == 'codon_start'):
                                    feature_codon_start = stri
                                elif (quali == 'translation'):
                                    if stri:
                                        # print type(stri), stri
                                        if len(stri[0]) > 0:
                                            seq_protein = stri[0]
                                            # print type(seq_protein), len(seq_protein), 'len protein', seq_protein
                                elif (quali == 'protein_id'):
                                    list_protein_ids = stri
                                elif (quali == 'lab_host'):
                                    feature_lab_host = stri
                                elif (quali == 'clone'):
                                    feature_clone = stri
                                elif (quali == 'note'):
                                    feature_note = stri
                                elif (quali == 'db_xref'):
                                    feature_db_xref = stri
                                elif (quali == 'tissue_type'):
                                    feature_tissue_type = stri
                                elif (quali == 'clone_lib'):
                                    feature_clone_lib = stri

                                if (printDetails): 
                                    print (offset*k), count2, ')', quali, stri
                            
                        else:
                            if (printDetails): 
                                print '********** unknown feature', feature
                            
                    indiv =  [self.basic.tira_aspas(indv_name),
                              self.basic.tira_aspas(indv_id),
                              self.basic.tira_aspas(indv_description),
                              self.basic.tira_aspas( str(version) ),
                              self.basic.tira_aspas(str_source),
                              self.basic.tira_aspas(str(list_taxonomy)),
                              self.basic.tira_aspas(str(list_keywords)),
                              self.basic.tira_aspas(str_date),
                              self.basic.tira_aspas(str_organism),
                              self.basic.tira_aspas(str_gi),
                              
                              self.basic.tira_aspas(str(list1_location)),
                              self.basic.tira_aspas(ref1_authors),
                              self.basic.tira_aspas(ref1_consortium),
                              self.basic.tira_aspas(ref1_journal),
                              self.basic.tira_aspas(ref1_medline),
                              self.basic.tira_aspas(ref1_pubmed),
                              self.basic.tira_aspas(ref1_comment),
                              
                              self.basic.tira_aspas(str(list2_location)),
                              self.basic.tira_aspas(ref2_authors),
                              self.basic.tira_aspas(ref2_consortium), 
                              self.basic.tira_aspas(ref2_journal),
                              self.basic.tira_aspas(ref2_medline),
                              self.basic.tira_aspas(ref2_pubmed),
                              self.basic.tira_aspas(ref2_comment),
                              
                              self.basic.tira_aspas(feature_strain), 
                              self.basic.tira_aspas(feature_mol_type),
                              self.basic.tira_aspas(feature_organism),
                              self.basic.tira_aspas(feature_db_xref),
                              self.basic.tira_aspas(list_genes),
                              self.basic.tira_aspas(list_products),
                              feature_codon_start,
                              seq_protein,
                              self.basic.tira_aspas(list_protein_ids),
                              self.basic.tira_aspas(feature_lab_host),
                              self.basic.tira_aspas(feature_clone),
                              self.basic.tira_aspas(feature_note),
                              self.basic.tira_aspas(feature_db_xref),
                              self.basic.tira_aspas(feature_tissue_type),
                              self.basic.tira_aspas(feature_clone_lib) 
                              ]
            
                    if printResume:
                        print '\n---- Resume --------------------------'
                        print 'organism', str_organism

                        if list_protein_ids:
                            print 'Protein ID'
                            for i in range(len(list_protein_ids)):
                                print offset, i+1,')',list_protein_ids[i]

                        print 'indv_id', indv_id
                        print '---------------------------------------'

                        if list_products:
                            print 'Product' 
                            for i in range(len(list_products)):
                                print offset, i+1,')',list_products[i]

                        if list_genes:
                            print 'Gene' 
                            for i in range(len(list_genes)):
                                print offset, i+1,')',list_genes[i]
                        
                        print '---------------------------------------'
                        
                        if list_taxonomy:
                            print 'taxonomy' 
                            for i in range(len(list_taxonomy)):
                                print offset, i+1,')',list_taxonomy[i]
                            
                        if feature_db_xref:
                            print 'db_xref' 
                            for i in range(len(feature_db_xref)):
                                print offset, i+1,')',feature_db_xref[i]
                        
                        if feature_mol_type:
                            print 'mol_type' 
                            for i in range(len(feature_mol_type)):
                                print offset, i+1,')',feature_mol_type[i]
                        print '---------------------------------------'
                        print 'Protein seq', seq_protein
                        print '---------------------------------------'
                        print 'DNA seq', seq_dna
                        print '----------------- end ------------------\n'
            except:
                print 'sequence error indv_id: ', indv_id
               
            handle.close()
            print 'Aquivo fechado.'
        
        print "\n--------------------------------------------------------"
        print 'End parsing.'


class NCBI_Read_Write:
    def __init__(self, db):
        self.db  = db
        self.acc = None
        self.genId = None
        self.specId = None
        
        self.basic = bioClass.Basic()
        
        Entrez.email = "flalix@gmail.com"
        print 'NCBI_Read_Write instantiated.'

    def openConnection(self):
        try:
            self.acc = cCrud.Crud()
            
            if (not self.acc.crud_connect()):
                print 'Error connection Wallace database1'
                return False
                
        except:
            print 'Error connection Wallace database2'
            return False

        print 'Database connected.'
        return True

    def queryWallace(self, spec_name, gene_searched, gene, dna_prot='PROT'):
        if not self.acc:
            print 'There is no connection.'
            return []  
            
        ret, rows = self.selectSequence(spec_name, gene_searched, gene, dna_prot)
              
        if not ret:
            print 'Tere are no rows'
            return []
        
        print 'read', len(rows), 'rows'
        return rows 
    
    def selectSequence(self, spec_name, gene_searched, gene, dna_prot = 'PROT'):
        if spec_name:
            sqlWhere = "upper(spec.name) = '" + spec_name.upper() + "'"
        else: 
            sqlWhere = None
                
        if gene_searched:
            if not sqlWhere:
                sqlWhere = "   upper(gene_searched) = '" + gene_searched.upper() + "' "
            else:
                sqlWhere = sqlWhere + "   and upper(gene_searched) = '" + gene_searched.upper() + "' "
        
        if gene:
            if not sqlWhere:
                sqlWhere = "   upper(feature_gene) = '" + gene.upper() + "'"
            else:
                sqlWhere = sqlWhere + "   and upper(feature_gene) = '" + gene.upper() + "'"
        
        if dna_prot == "DNA":
            query =  \
            "select spec.cd, spec.name, gene_searched, " + \
            "       feature_gene, indv_id, sequence, " + \
            "       feature_db_xref db_xref, indv_name, indv_description " + \
            "from  species_individuals sindv " + \
            "join  species spec on spec.id = sindv.species " + \
            "where " + sqlWhere
                    
        else:
            query = \
            "select spec.cd, spec.name, gene_searched, " + \
            "       feature_gene, indv_id, feature_translation sequence, " + \
            "       feature_db_xref db_xref, indv_name, indv_description " + \
            "from  species_individuals sindv " + \
            "join  species spec on spec.id = sindv.species " + \
            "where " + sqlWhere
            
        print '**', query

        return self.acc.crud_query(query)
    
    
    
    def defineDrosophilaGenus(self):
        self.domainID = self.acc.createDomain('EUKARYA', 'Eukarya')
        if self.domainID <= 0:
            print 'Domain not defined'
            return False
        print 'domainID', self.domainID
        
        self.kingId = self.acc.createKingdom(self.domainID, 'ANIMALIA', 'Animalia')
        if self.kingId <= 0:
            print 'Kingdom not defined'
            return False
        print 'kingId', self.kingId

        self.philId = self.acc.createPhylum(self.kingId, 'ARTHOPODA', 'Arthopoda')
        if self.philId <= 0:
            print 'Philo not defined'
            return False
        print 'philId', self.philId

        self.classId = self.acc.createClass(self.philId, 'HEXAPODA', 'Hexapoda')
        if self.classId <= 0:
            print 'Class not defined'
            return False
        print 'classId', self.classId

        self.ordId = self.acc.createOrder(self.classId, 'DIPTERA', 'Diptera')
        if self.ordId <= 0:
            print 'Order not defined'
            return False
        print 'ordId', self.ordId

        self.famId = self.acc.createFamily(self.ordId, 'DROSOPHILA', 'Drosophilidae')
        if self.famId <= 0:
            print 'Family not defined'
            return False
        print 'famId', self.famId

        self.genId = self.acc.createGenus(self.famId, 'DROSOPHILA', 'Drosophila')
        if self.genId <= 0:
            print 'Genus not defined'
            return False
        print 'genId', self.genId

        return True
                
    def writeDB(self, gene, idList, printResume=False, printDetails=False, showComment=False):
        
        if not self.genId:
            if not self.defineDrosophilaGenus():
                print "Genus Drosophila is not defined"
                return False 

        if not gene:
            print "Gene searched not defined1"
            return False 

        if len(gene) == 0:
            print "Gene searched not defined2"
            return False 

        lenIdList = len(idList)
        print lenIdList, 'species selected: ', idList 
        
        for i in range( lenIdList ):
            
            access_ID = idList[i]
            
            print "--------------------------------------------------------"
            print 'access_ID', access_ID
            
            try:
                handle = Entrez.efetch(db=self.db, rettype="gb", id=access_ID, retmode="text")
                
            except:
                print '   !!!!! error in access_ID   !!!!! Stop parsing !!!!', access_ID
                return False

            offset = '  '
            
            try:
                # print 'Parsing NCBI seq ...'
                # print SeqIO.parse(handle, "gb")
                
                # for seq_record in Entrez.parse(handle):
                for seq_record in SeqIO.parse(handle, "gb"):
                    # print 'seq_record....'
                    print "### Parsing '%s', length %i, with %i features" % (seq_record.name, len(seq_record), len(seq_record.features))
        
                    indv_description = seq_record.description
                    indv_name = seq_record.name
                    indv_id = seq_record.id
                    seq_dna = seq_record.seq
                    # print len(seq_dna), 'len dna', seq_dna

                    version = 0
                    str_source = ''
                    comment = ''

                    list_taxonomy = []
                    list_keywords = []
                    list_accessions = []

                    str_date_file_division = ''
                    str_date = ''
                    str_organism = ''
                    str_gi = ''

                    feature_strain = ''
                    feature_mol_type = ''
                    feature_organism = ''
                    feature_db_xref = ''
                    list_genes = []
                    list_products = []
                    feature_codon_start = ''
                    seq_protein = ''   # is a list !!
                    list_protein_ids = ''
                    feature_lab_host = ''
                    
                    feature_clone = ''
                    feature_note = ''
                    feature_db_xref = ''
                    feature_tissue_type = ''
                    feature_clone_lib = ''
                            
                    list1_location = []
                    ref1_authors = None
                    ref1_consortium = None
                    ref1_journal = None
                    ref1_medline = None
                    ref1_pubmed = None
                    ref1_comment = None
                    
                    list2_location = []
                    ref2_authors = None
                    ref2_consortium  = None
                    ref2_journal = None
                    ref2_medline = None
                    ref2_pubmed = None
                    ref2_comment = None   
                               
                              
                    i = 1
                    if (printDetails):
                        print 'Keys', seq_record.annotations.keys()
                        
                    for annot in seq_record.annotations:
            
                        j = 1
                        if (annot == 'references'):
                            if (printDetails):
                                print '\n=== Refereces ========================'
                                            
                            annotNumber = 1
                            refObjects = seq_record.annotations[annot]
            
                            for ref in refObjects:
                                # print ref
                                
                                if (annotNumber == 1):
                                        list1_location = ref.location
                                        ref1_authors = ref.authors
                                        ref1_consortium = ref.consrtm
                                        ref1_journal = ref.journal
                                        ref1_medline = ref.medline_id
                                        ref1_pubmed = ref.pubmed_id
                                        ref1_comment = ref.comment
                                        
                                        if (printDetails):
                                            print (offset*(i+1)), str(j) + '. reference1:'
                                            if list1_location:
                                                print (offset*(i+2)), 'location:', list1_location
                                            if ref1_authors:
                                                print (offset*(i+2)), 'authors:', ref1_authors
                                            if ref1_consortium:
                                                print (offset*(i+2)), 'consortium:', ref1_consortium
                                            if ref1_journal:
                                                print (offset*(i+2)), 'journal:', ref1_journal
                                            if ref1_medline:
                                                print (offset*(i+2)), 'medline:', ref1_medline
                                            if ref1_pubmed:
                                                print (offset*(i+2)), 'pubmed:', ref1_pubmed
                                            if ref1_comment:
                                                print (offset*(i+2)), 'comment:', ref1_comment 
                                            
                                            print '------- end ref1 ---------------\n'
                                                                    
                                elif (annotNumber == 2):        
                                        list2_location = ref.location
                                        ref2_authors = ref.authors
                                        ref2_consortium = ref.consrtm
                                        ref2_journal = ref.journal
                                        ref2_medline = ref.medline_id
                                        ref2_pubmed = ref.pubmed_id
                                        ref2_comment = ref.comment
                                        
                                        if (printDetails):
                                            print (offset*(i+1)), str(j) + '. reference2:'
                                            if list2_location:
                                                print (offset*(i+2)), 'location:', list2_location
                                            if ref2_authors:
                                                print (offset*(i+2)), 'authors:', ref2_authors
                                            if ref2_consortium:
                                                print (offset*(i+2)), 'consortium:', ref2_consortium
                                            if ref2_journal:
                                                print (offset*(i+2)), 'journal:', ref2_journal
                                            if ref2_medline:
                                                print (offset*(i+2)), 'medline:', ref2_medline
                                            if ref2_pubmed:
                                                print (offset*(i+2)), 'pubmed:', ref2_pubmed
                                            if ref2_comment:
                                                print (offset*(i+2)), 'comment:', ref2_comment 

                                            print '------- end ref2 ---------------\n'
                                        
                                annotNumber += 1
                                j += 1
                        else:
                            if (printDetails):
                                print '\n==== Annotations ======================='
                            
                            varAny = seq_record.annotations[annot]
                            
                            if (annot == 'comment'):
                                comment = varAny
                                if showComment:
                                    print 'comment ', version
                                
                            elif (annot == 'sequence_version'):
                                version = varAny
                                if (printDetails):
                                    print 'sequence_version ', version, type(version)
                            elif (annot == 'source'):
                                str_source = varAny
                                if (printDetails):
                                    print offset, 'source', varAny
                            elif (annot == 'taxonomy'):
                                list_taxonomy = varAny
                                if (printDetails):
                                    print offset, 'taxonomy', varAny
                            elif (annot == 'date'):
                                str_date = varAny
                                if (printDetails):
                                    print offset, 'date', varAny
                            elif (annot == 'organism'):
                                str_organism = varAny
                                if (printDetails):
                                    print offset, 'organism', varAny
                            elif (annot == 'gi'):
                                str_gi = varAny
                                if (printDetails):
                                    print offset, 'gi', varAny
                            elif (annot == 'accessions'):
                                list_accessions = varAny
                                if (printDetails):
                                    print offset, 'accessions', varAny
                            elif (annot == 'data_file_division'):
                                str_date_file_division = varAny
                                if (printDetails):
                                    print offset, 'data_file_division', varAny
                            else:
                                if (printDetails):
                                    print offset, annot, varAny
                                
                    if (printDetails):    
                        print '\n----------------------------'
                        print "%i features" % len(seq_record.features)
                        # print "%i sub-features," % len(seq_record.SubFeatures)
                      
            
                    for feature in seq_record.features:
                        if (printDetails): 
                            print '\n', offset, str(i) + '. feature:'
                        i += 1
                        
                        j = 2
                        count = 1
                        if (feature.type):
                            stri = feature.qualifiers.get(feature.type)
                            if (stri != None):
                                if (printDetails): 
                                    print (offset*j), count,') type', feature.type, stri
                                count += 1
            
                        if (feature.location):
                            stri = feature.qualifiers.get(feature.location)
                            if (stri != None):
                                if (printDetails): 
                                    print (offset*j), count,') location',  stri
                                count += 1
            
                        if (feature.strand):
                            stri = feature.qualifiers.get(feature.strand)
                            if (stri != None):
                                if (printDetails): 
                                    print (offset*j), count,') strand', feature.strand, stri
                                count += 1
            
                        if (feature.qualifiers):
                            if (printDetails): 
                                print (offset*j), count,') qualifiers:'

                            count += 1
            
                            k = 3
                            count2 = 0
                            for quali in feature.qualifiers:
                                stri = feature.qualifiers.get(quali)
                                count2 += 1
                               
                                if (quali == 'strain'):
                                    feature_strain = stri
                                elif (quali == 'mol_type'):
                                    feature_mol_type = stri
                                elif (quali == 'organism'):
                                    feature_organism = stri
                                elif (quali == 'db_xref'):
                                    feature_db_xref = stri
                                elif (quali == 'gene'):
                                    list_genes = stri
                                elif (quali == 'product'):
                                    list_products = stri
                                elif (quali == 'codon_start'):
                                    feature_codon_start = stri
                                elif (quali == 'translation'):
                                    if stri:
                                        # print type(stri), stri
                                        if len(stri[0]) > 0:
                                            seq_protein = stri[0]
                                            # print type(seq_protein), len(seq_protein), 'len protein', seq_protein
                                elif (quali == 'protein_id'):
                                    list_protein_ids = stri
                                elif (quali == 'lab_host'):
                                    feature_lab_host = stri
                                elif (quali == 'clone'):
                                    feature_clone = stri
                                elif (quali == 'note'):
                                    feature_note = stri
                                elif (quali == 'db_xref'):
                                    feature_db_xref = stri
                                elif (quali == 'tissue_type'):
                                    feature_tissue_type = stri
                                elif (quali == 'clone_lib'):
                                    feature_clone_lib = stri

                                if (printDetails): 
                                    print (offset*k), count2, ')', quali, stri
                            
                        else:
                            if (printDetails): 
                                print '********** unknown feature', feature
                            
                    indiv =  [gene.upper().strip(),
                              self.basic.tira_aspas(indv_name),
                              self.basic.tira_aspas(indv_id),
                              self.basic.tira_aspas(indv_description),
                              self.basic.tira_aspas( str(version) ),
                              self.basic.tira_aspas(str_source),
                              self.basic.tira_aspas(str(list_taxonomy)),
                              self.basic.tira_aspas(str(list_keywords)),
                              self.basic.tira_aspas(str_date),
                              self.basic.tira_aspas(str_organism),
                              self.basic.tira_aspas(str_gi),
                              
                              self.basic.tira_aspas(str(list1_location)),
                              self.basic.tira_aspas(ref1_authors),
                              self.basic.tira_aspas(ref1_consortium),
                              self.basic.tira_aspas(ref1_journal),
                              self.basic.tira_aspas(ref1_medline),
                              self.basic.tira_aspas(ref1_pubmed),
                              self.basic.tira_aspas(ref1_comment),
                              
                              self.basic.tira_aspas(str(list2_location)),
                              self.basic.tira_aspas(ref2_authors),
                              self.basic.tira_aspas(ref2_consortium), 
                              self.basic.tira_aspas(ref2_journal),
                              self.basic.tira_aspas(ref2_medline),
                              self.basic.tira_aspas(ref2_pubmed),
                              self.basic.tira_aspas(ref2_comment),
                              
                              self.basic.tira_aspas(feature_strain), 
                              self.basic.tira_aspas(feature_mol_type),
                              self.basic.tira_aspas(feature_organism),
                              self.basic.tira_aspas(feature_db_xref),
                              self.basic.tira_aspas(list_genes),
                              self.basic.tira_aspas(list_products),
                              self.basic.tira_aspas(feature_codon_start),
                              self.basic.tira_aspas(seq_protein),
                              self.basic.tira_aspas(seq_dna),
                              self.basic.tira_aspas(list_protein_ids) 
                              # self.basic.tira_aspas(feature_lab_host),
                              # self.basic.tira_aspas(feature_clone),
                              # self.basic.tira_aspas(feature_note),
                              # self.basic.tira_aspas(feature_tissue_type),
                              # self.basic.tira_aspas(feature_clone_lib) 
                              ]
            
                    if printResume:
                        print '\n---- Resume --------------------------'
                        print 'organism', str_organism

                        if list_protein_ids:
                            print 'Protein ID'
                            for i in range(len(list_protein_ids)):
                                print offset, i+1,')',list_protein_ids[i]

                        print 'indv_id', indv_id
                        print '---------------------------------------'

                        if list_products:
                            print 'Product' 
                            for i in range(len(list_products)):
                                print offset, i+1,')',list_products[i]

                        if list_genes:
                            print 'Gene' 
                            for i in range(len(list_genes)):
                                print offset, i+1,')',list_genes[i]
                        
                        print '---------------------------------------'
                        
                        if list_taxonomy:
                            print 'taxonomy' 
                            for i in range(len(list_taxonomy)):
                                print offset, i+1,')',list_taxonomy[i]
                            
                        if feature_db_xref:
                            print 'db_xref' 
                            for i in range(len(feature_db_xref)):
                                print offset, i+1,')',feature_db_xref[i]
                        
                        if feature_mol_type:
                            print 'mol_type' 
                            for i in range(len(feature_mol_type)):
                                print offset, i+1,')',feature_mol_type[i]
                        print '---------------------------------------'
                        print 'Protein seq', seq_protein
                        print '---------------------------------------'
                        print 'DNA seq', seq_dna
                    else:
                        print '\n---- No Resume --------------------------'
                        print 'organism', str_organism
                        print 'indv_id', indv_id
                        if list_genes:
                            print 'Gene', list_genes[0]
                        else:
                            print 'Gene None?'

                        
                    print '----------------- end ------------------\n'

                    try:
                        specDs = str_organism.replace('Drosophila','').strip()
                        specCD = specDs.upper()[0:10]
                        
                        self.specId = self.acc.createSpecies(self.genId, specCD, specDs)
                        if self.specId <= 0:
                            print 'Species ' + specDs + ' is not defined'
                            handle.close()
                            self.acc.crud_close() 
                            return False
                        
                        # print 'specId', self.specId
                        
                        indivDbID = self.acc.createIndivSpecies(self.specId, indv_id, indiv, showMessage=False)
                        print 'wrote indv_id', indv_id, ' db ID:', indivDbID
                          
                    except:
                        print 'Error saving indv_id:', indv_id, specDs
                        handle.close()
                        self.acc.crud_close()
                        return False
                                                        
            except:
                print 'sequence error indv_id: ', indv_id
                handle.close()
                self.acc.crud_close()
                return False
                    
            handle.close()
            print 'Aquivo fechado.'
        
        print "\n--------------------------------------------------------"
        self.acc.crud_close()
        print 'End saving Drosophilas''s indivs.'
        return True
        

    
    
        
                            
                    



    