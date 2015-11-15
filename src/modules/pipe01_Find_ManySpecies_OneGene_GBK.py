'''
Created on 02/09/2013
Update  on 01/07/2014
Update  on 05/10/2015 - get_params()

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica
'''
from Bio import Entrez
from Bio import SeqIO
import os

'''
  for adh search: Capsella , Drosophila , Calystegia, Staphylococcus
'''

class Pipe():
    def __init__(self, desk):
        
        self.failure = True
        self.error_msg = ''
        self.desk = desk
        

        try:
            desk.get_params()
            desk.showGraph = False
            
            if desk.geneList_var.get() == '':
                if self.desk.gene == '':
                    self.desk.geneList = ['']
                else:
                    self.desk.geneList = [self.desk.gene]
                
                
            if desk.specieList_var.get() == '':
                specieList = ['']
            else:    
                specieList = desk.specieList_var.get().split(',')


            self.filenameGBK = desk.rootFasta+desk.output_filename_gbk

            ''' obligatory '''
            Entrez.email = desk.email_var.get()
            
        except:
            self.error_msg = 'Could not get parameters.'
            return

        if 'xxx' in Entrez.email:
            'Please define your email'
            self.error_msg = 'Please define your email'
            return
        

        if self.desk.organism == '':
            self.error_msg = 'Define the organism.'
            return
        
        if not self.desk.geneList and self.desk.gene == '' and self.desk.title == '':
            self.error_msg = 'Define the default gene or title.'
            return


        logDic = {}
        filename = desk.rootFasta + desk.logGbkFilename_var.get()
        
        if os.path.isfile(filename):
            try:
                f = open(filename, 'r')
                
                print 'reading %s'%(filename)
                for line in f:
                    idNCBI,_ = line.split('\n')
                    logDic[idNCBI] = idNCBI
                    
            finally:
                f.close
        
        
        print self.filenameGBK
        
        for this_gene in self.desk.geneList:
            for species in specieList:

                if self.desk.organism.lower() == 'any':
                    searchTerm = ''
                else:
                    searchTerm = '(%s [organism])' % (self.desk.organism)
                
                if this_gene <> '':
                    if this_gene.upper() != 'ANY':
                        if searchTerm == '':
                            searchTerm += '(%s[Gene])' % (this_gene)
                        else:
                            searchTerm += ' and (%s[Gene])' % (this_gene)

                if self.desk.title != '': 
                    searchTerm += self.desk.title + '[Title]'
                    
                '''
                if title != '': 
                    titleList = title.split(' ')
                    title = ''
                    
                    and_or = False
                    
                    for w in titleList:
                        w = w.strip()
                        if title == '':
                            title += w # +  '[Title]'
                        else:
                            if w.lower() == 'and':
                                title += ' and '
                                and_or = True
                            elif w.lower() == 'or':
                                title += ' or '
                                and_or = True
                            else:
                                if and_or:
                                    title += w # +  '[Title]'
                                    and_or = False
                                else:
                                    title += ' and ' + w # + '[Title]'
                    
                                        
                    if searchTerm == '':
                        searchTerm += ' (%s)' % (title)
                    else:
                        searchTerm += ' and (%s)' % (title)
                '''
                if self.desk.completeSeq == 'no':
                    searchTerm += ' and not (complete sequence[Title])'
                elif self.desk.completeSeq == 'yes':
                    searchTerm += ' and (complete sequence[Title])'
                    
                if self.desk.completeGen == 'no':
                    searchTerm += ' and not (complete genome[Title])'
                elif self.desk.completeGen == 'yes':
                    searchTerm += ' and (complete genome[Title])'
                    
                if self.desk.chromosome == 'no':
                    searchTerm += ' and not (chromosome[Title] or chromosomal[Title])'
                elif self.desk.chromosome == 'yes':
                    searchTerm += ' and (chromosome[Title])'
                    
                if self.desk.shotGun == 'no':
                    searchTerm += ' and not (shotgun[Title])'
                elif self.desk.shotGun == 'yes':
                    searchTerm += ' and (shotgun[Title])'
                    
                if self.desk.contig == 'no':
                    searchTerm += ' and not (genomic contig[Title])'
                elif self.desk.contig == 'yes':
                    searchTerm += 'and (genomic contig[Title] '
                    
                if self.desk.superContig == 'no':
                    searchTerm += ' and not (supercontig [Title])'
                elif self.desk.superContig == 'yes':
                    searchTerm += ' and (supercontig [Title])'
                    
                if self.desk.partialCDS == 'no':
                    searchTerm += ' and not (partial cds [Title])'
                elif self.desk.partialCDS == 'yes':
                    searchTerm += ' and (partial cds [Title])'
        
                if self.desk.cds == 'no':
                    searchTerm += ' and not (cds [Title])'
                elif self.desk.cds == 'yes':
                    searchTerm += ' and (cds [Title])'

                if self.desk.mrna == 'no':
                    searchTerm += ' and not (mRNA [Title])'
                elif self.desk.mrna == 'yes':
                    searchTerm += ' and (mRNA [Title])'
                    
                searchTerm += ' and not (synthetic [Title])'

                print '\n----------------------------'
                desk.showmsg_obs("Search Term: %s"%(searchTerm))
                
                handle = Entrez.esearch(db=self.desk.db, term=searchTerm, rettype='gb', retmax=self.desk.retmax)
                record = Entrez.read(handle)
                lista = record["IdList"]
   
                j = 0
                for idNCBI in lista:
                    j += 1
                    records = None
                    seq_record = None
                    
                    if idNCBI == "20279117":
                        pass
                    
                    if species == '':
                        stri = "%s %i / %i, gene %s "%(idNCBI, j, len(lista), this_gene)
                    else:
                        stri = "%s %i / %i, species %s  gene %s "%(idNCBI, j, len(lista), species, this_gene)
                        
                    geneij = this_gene + ' - ' + idNCBI
                    
                    if idNCBI in logDic.keys():
                        stri += '### This seq is already saved'
                        desk.showmsg_obs(stri)
                        continue
        
        
                    try:
                        handle = Entrez.efetch(db=self.desk.db, rettype="gb", retmode="text", id=idNCBI)
                        ''' work as long as the file contains at least one record '''
                        records = None
                        seq_record = SeqIO.read(handle, "gb")
                    except:
                        '''
                        http://biopython.org/DIST/docs/api/Bio.SeqIO-module.html
                        For non-interlaced files (e.g. Fasta, GenBank, EMBL) with multiple records using a 
                        sequence iterator can save you a lot of memory (RAM). There is less benefit for interlaced file formats (e.g. most multiple alignment file formats). However, an iterator only lets you access the records one by one.
                        If you want random access to the records by number, turn this into a list:
                        '''
                        try:
                            records = list(SeqIO.parse(handle, "gb"))     
                        except:
                            desk.showmsg_obs("Reading error, problems with NCBI connection. Try again")
                            return
                    finally:
                        handle.close()
                    
                    its_ok = False
                    if records:
                        desk.showmsg_obs(stri + '>>> found records')
                        for seq_record in records:
                            try:
                                if not os.path.isfile(self.filenameGBK):
                                    handle = open(self.filenameGBK,"w")
                                else:
                                    handle = open(self.filenameGBK,"a")
            
                                SeqIO.write(seq_record, handle, "genbank")
                                its_ok = True
                            except:
                                its_ok = False
                            finally:
                                handle.close()
        
                    else:
                        desk.showmsg_obs(stri)
                        if seq_record:
                            try:
                                handle = open(self.filenameGBK,"a")
                                SeqIO.write(seq_record, handle, "genbank")
                                its_ok = True
                            except:
                                print 'Could save record %s %s'%(seq_record.id, geneij)
                                its_ok = False
                            finally:
                                handle.close()
                                                  
                        else:
                            print '  ****  ID:', idNCBI, 'not found in NCBI.'
            
                        if its_ok:
                            try:
                                f = open(filename, 'a')
                                f.write(idNCBI + '\n')
                            except:
                                print 'Could not write %s %s %s'%(species, this_gene, idNCBI)
                                
                            finally:
                                f.close()  
        
        self.failure = False
        return
