'''
Created on 04/02/2012
Last Edited on 05/02/2012

@author: Flavio Lichtenstein
'''

import MySQLdb as mdb
con = None

''' 
mysql:  root:admin localhost:3307 centro
http://www.lfd.uci.edu/~gohlke/pythonlibs/#mysql-python
http://mysql-python.sourceforge.net/MySQLdb-1.2.2/private/_mysql-module.html
'''

class Crud:
    def __init__(self):
        self.con = None
        self.database = 'walace'
        self.user = 'root'
        self.pwd  = 'admin'
        
        
    def crud_connect(self, showMessage=True):
        if (self.con != None):
            return True
        
        try:
            self.con = mdb.connect(host="localhost", \
                    user=self.user, passwd=self.pwd, db=self.database)
            if showMessage:
                print 'MySql Connection stabilished.'
            return True
        
        except mdb.Error, e:
            print "Error %d: %s" % (e.args[0], e.args[1])
            return False
        
            
    def crud_query(self, query):   
        try:
            query = query.lstrip().rstrip().upper()

            # print 'crud_query', query[0:6], "'" + query + "'"
            cur = self.con.cursor()
            cur.execute(query)
        
            rows = []
            if (query[0:6] == 'SELECT'):
                rows = cur.fetchall()
                
                # if (len(rows) == 0):
                #     print 'there are no lines.'
                # may be into cursor
                cur.close
                # print 'query', query
                # print '*** rows', rows
                if len(rows) == 0:
                    ret = False
                else:
                    ret = True
            else:
                self.con.commit
                cur.close
                ret = True
            
        except mdb.Error, e:
            cur.close
            print "Error %d: %s" % (e.args[0], e.args[1])
            ret = False
            
        return ret, rows

    
    def crud_insert(self, table, fields, values, showmessage=False):
       
        try:
            with self.con: 
                if (table == None):
                    print 'Error: Name from the table, please.'
                    return False
                if (values == None):
                    print 'Error: values does not exist, please.'
                    return False
                
                cur = self.con.cursor()

                if (fields == None):
                    query = "INSERT INTO " + table + " VALUES (" + values + ")"
                else:
                    query = "INSERT INTO " + table + " (" + fields + ") VALUES (" + values + ")"

                if showmessage:
                    print '*** insert: ', query
                
                try:
                    cur.execute(query)
                except:
                    print 'Error: exception'
                    return False
            
                return True
            
        except mdb.Error, e:
          
            print "Error %d: %s" % (e.args[0], e.args[1])
            return False 


    def crud_exists(self, table, field, value):
        try:
            with self.con: 
                if (table == None):
                    print 'Error: Name from the table, please.'
                    return -1
                if (field == None):
                    print 'Error: field does not exist, please.'
                    return -1
                if (value == None):
                    print 'Error: value does not exist, please.'
                    return -1
                
                query = "SELECT ID FROM " + table + " WHERE " + field + " = '" + value + "'"
                # print 'query -->: ', query
                ret, row = self.crud_query(query)
                
                if not ret:
                    return -1

                if row == None:
                    return -1
                
                if not row[0]:
                    return -1
                
                reg = row[0]
                # print 'ok crud_exists', reg[0], row, 'query -->: ', query

                return reg[0] 
            
        except mdb.Error, e:
          
            print "Error %d: %s" % (e.args[0], e.args[1])
            return -1
                

    def crud_exists_father(self, table, fatField, fatValue, field, value, showmessage=False):
        try:
            with self.con: 
                if (table == None):
                    print 'Error: Name from the table, please.'
                    return -1
                if (field == None):
                    print 'Error: field does not exist, please.'
                    return -1
                if (value == None):
                    print 'Error: value does not exist, please.'
                    return -1
                if (fatField == None):
                    print 'Error: father Field does not exist, please.'
                    return -1
                if (fatValue == None):
                    print 'Error: father Value does not exist, please.'
                    return -1

                query = "SELECT ID FROM " + table + \
                        " WHERE " + fatField + " = "  + str(fatValue) + " AND " \
                                  + field    + " = '" + value + "'"

                if (showmessage):
                    print "crud_exists_father - query: ", query
                    
                ret, row = self.crud_query(query)

                if not ret:
                    return -1

                if row == None:
                    return -1
                
                if not row[0]:
                    return -1
                
                reg = row[0]
                # print 'ok crud_exists_father', reg[0], row, 'query -->: ', query

                return reg[0] 
 
            
        except mdb.Error, e:
          
            print "Error %d: %s" % (e.args[0], e.args[1])
            return -1
        
    def crud_delete(self, table, field, symbol, regId):

        try:
            if (table == None):
                print 'Error: Name from the table, please.'
                return -1
            if (field == None):
                print 'Error: field does not exist, please.'
                return -1
            if (regId == None):
                print 'Error: regId does not exist, please.'
                return -1
            if (regId <= 0):
                print 'Error: regId must be > 0.'
                return -1                
            
            query = "DELETE FROM " + table + " WHERE " + field + " " + symbol + ' ' + str(regId)
            # print query
            cur = self.con.cursor()
            cur.execute(query)  
            self.con.commit
            cur.close()

            
            return True
            
        except mdb.Error, e:
            self.con.rollback()
            print "Error %d: %s" % (e.args[0], e.args[1])
            return False
                
   
    def crud_close(self):    
        if self.con:
            self.con.close()
            
    
    def createDomain(self, cd, name):
        if (not self.crud_connect()):
            print 'nao foi possivel conectar em Walace'

        ret = self.crud_exists('DOMAIN', 'CD', cd)
        
        if (ret <= 0):
            self.crud_insert('domain', 'cd, name, created_by', "'" + cd + "', '" + name + "', 99")
            ret = self.crud_exists('DOMAIN', 'CD', cd)
            
        # print 'ret',ret   
        return ret
    

    def createKingdom(self, fatId, cd, name):
        if (not self.crud_connect()):
            print 'nao foi possivel conectar em Walace'
            return -1

        ret = self.crud_exists_father('KINGDOM', 'domain', fatId, 'CD', cd)
        
        if (ret <= 0):
            self.crud_insert('KINGDOM', 'domain, cd, name, created_by', \
                             str(fatId) + ", '" + cd + "', '" + name + "', 99")
            ret = self.crud_exists_father('KINGDOM', 'domain', fatId, 'CD', cd)
            
        # print 'ret',ret   
        return ret
    
    
    def createPhylum(self, fatId, cd, name):
        if (not self.crud_connect()):
            print 'nao foi possivel conectar em Walace'
            return -1

        ret = self.crud_exists_father('PHYLUM', 'kingdom', fatId, 'CD', cd)
        
        if (ret <= 0):
            self.crud_insert('PHYLUM', 'kingdom, cd, name, created_by', \
                             str(fatId) + ", '" + cd + "', '" + name + "', 99")
            ret = self.crud_exists_father('PHYLUM', 'kingdom', fatId, 'CD', cd)
            
        # print 'ret',ret   
        return ret


    def createClass(self, fatId, cd, name):
        if (not self.crud_connect()):
            print 'nao foi possivel conectar em Walace'
            return -1
     
        ret = self.crud_exists_father('phylum_class', 'phylum', fatId, 'CD', cd)
        
        if (ret <= 0):
            print 'class --> phylum', fatId
            self.crud_insert('phylum_class', 'phylum, cd, name, created_by', \
                             str(fatId) + ", '" + cd + "', '" + name + "', 99")
            ret = self.crud_exists_father('phylum_class', 'phylum', fatId, 'CD', cd)
            
        # print 'ret',ret   
        return ret

    
    def createOrder(self, fatId, cd, name):
        if (not self.crud_connect()):
            print 'nao foi possivel conectar em Walace'
            return -1
     
        ret = self.crud_exists_father('class_order', 'phylum_class', fatId, 'CD', cd)
        
        if (ret <= 0):
            print 'order --> class', fatId
            self.crud_insert('class_order', 'phylum_class, cd, name, created_by', \
                             str(fatId) + ", '" + cd + "', '" + name + "', 99")
            ret = self.crud_exists_father('class_order', 'phylum_class', fatId, 'CD', cd)
            
        # print 'ret',ret   
        return ret
    
    
    def createFamily(self, fatId, cd, name):
        if (not self.crud_connect()):
            print 'nao foi possivel conectar em Walace'
            return -1
     
        ret = self.crud_exists_father('family', 'class_order', fatId, 'CD', cd)
        
        if (ret <= 0):
            print 'family --> order', fatId
            self.crud_insert('family', 'class_order, cd, name, created_by', \
                             str(fatId) + ", '" + cd + "', '" + name + "', 99")
            ret = self.crud_exists_father('family', 'class_order', fatId, 'CD', cd)
            
        # print 'ret',ret   
        return ret
    
    
       
    def createGenus(self, fatId, cd, name):
        if (not self.crud_connect()):
            print 'nao foi possivel conectar em Walace'
            return -1
     
        ret = self.crud_exists_father('genus', 'family', fatId, 'CD', cd)
        
        if (ret <= 0):
            print 'genus --> family', fatId            
            self.crud_insert('genus', 'family, cd, name, created_by', \
                             str(fatId) + ", '" + cd + "', '" + name + "', 99")
            ret = self.crud_exists_father('genus', 'family', fatId, 'CD', cd)
            
        # print 'ret',ret   
        return ret
    

    def createSpecies(self, fatId, cd, name):
        if (not self.crud_connect()):
            print 'nao foi possivel conectar em Walace'
            return -1
     
        ret = self.crud_exists_father('species', 'genus', fatId, 'CD', cd)
        
        if (ret <= 0):
            print 'species --> genus', fatId            
            self.crud_insert('species', 'genus, cd, name, created_by', \
                             str(fatId) + ", '" + cd + "', '" + name + "', 99")
            ret = self.crud_exists_father('species', 'genus', fatId, 'CD', cd)
            
        # print 'ret',ret   
        return ret
    

    def createIndivSpecies(self, specId, indv_id, indiv, showMessage=False):
        if (not self.crud_connect()):
            print 'nao foi possivel conectar em Walace'
            return -1
            
        ret = self.crud_exists_father('Species_Individuals', 'species', specId, 'indv_id', indv_id)
        
        if (ret > 0):
            print 'indiv exists: species', specId, 'indv_id', indv_id
            return ret
        
        fields = 'species, gene_searched, indv_name, indv_id, indv_description, indv_sequence_version,\
                  indv_source, indv_taxonomy, indv_keywords, indv_date, \
                  indv_organism, indv_gi, ref1_location, ref1_authors, \
                  ref1_consortium, ref1_journal, ref1_medline, ref1_pubmed, \
                  ref1_comment, ref2_location, ref2_authors, ref2_consortium, \
                  ref2_journal, ref2_medline, ref2_pubmed, ref2_comment, \
                  feature_strain, feature_mol_type, feature_organism, \
                  feature_db_xref, feature_gene, feature_product, \
                  feature_codon_start, feature_translation, sequence, feature_protein_id'
                          
                      
        values = str(indiv)
        
        values = values[1:len(values)-1]
        values = values.replace("[",'')
        values = values.replace("]",'')
        values = values.replace("'Null'",'Null')
        values = values.replace("''",'Null')

        # including geneus
        values = str(specId) + ', ' +  values
        
        if showMessage:
            print '\n--------------------------'
            print 'fields:  ', fields
            print '--------------------------'
            print 'values:  ', values
            print '--------------------------'
        

        if not self.crud_insert('Species_Individuals', fields, values):
            print 'could not save: species', specId, 'indv_id', indv_id
            return -1
        
        ret = self.crud_exists_father('Species_Individuals', 'species', specId, 'indv_id', indv_id)
            
        if showMessage:
            print 'saving: species', specId, 'indv_id', indv_id, 'reg = ', ret
        return ret
    

