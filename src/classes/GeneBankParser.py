'''
Created on 15/08/2013
Created on 26/04/2014 - %s gene
Updated on 19/03/2015 - params = dic_species_params[species] 

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica

http://biopython.org/wiki/SeqRecord
http://biopython.org/wiki/SeqIO
http://biopython.org/DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html
http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/
http://pythonadventures.wordpress.com/archives/ 
'''

from Bio import SeqIO
from Bio.Seq import Seq
# from Bio.Alphabet import IUPAC
from BioPythonClass import MySequence
from Bio.SeqRecord import SeqRecord
import os


class GeneBankParser:
    def __init__(self, desk):
        
        self.species = 'unknown'
        self.organism = desk.organism
        self.species = desk.species
        self.title = desk.title
        self.gene = desk.gene
        self.gene_title = desk.gene_title
        self.cutoffLength = desk.cutoffLength
        self.rootFasta = desk.rootFasta
        
        '''                      Drosophila_Adh '''
        self.filenameGBK = desk.rootFasta + desk.output_filename_gbk
    
        self.filename_GeneFasta = desk.rootFasta+'%s_all_Gene_%s_cutoff%i.fasta' % (desk.organism, desk.gene_title, desk.cutoffLength)
        self.filename_CDSFasta = desk.rootFasta+'%s_all_CDS_%s_cutoff%i.fasta' % (desk.organism, desk.gene_title, desk.cutoffLength)
        self.filename_ExonFasta = desk.rootFasta+'%s_all_Exon_%s_cutoff%i.fasta' % (desk.organism, desk.gene_title, desk.cutoffLength)
        self.filename_IntronFasta = desk.rootFasta+'%s_all_Intron_%s_cutoff%i.fasta' % (desk.organism, desk.gene_title, desk.cutoffLength)
        self.filename_ProtFasta = desk.rootFasta+'%s_all_Prot_%s_cutoff%i.fasta' % (desk.organism, desk.gene_title, desk.cutoffLength)
        
        self.qualifiers = ('Sourcself.cutoffLengthGene','mRNA','CDS','exon','intron',"5'UTR",'CAAT_signal','TATA_signal','promoter','variation','polyA_signal','polyA_site','misc_feature','repeat_region','precursor_RNA','gap','prim_transcript','mobile_element','LTR','3UTR','misc_difference','enhancer','mat_peptide','rRNA','misc_RNA','primer_bind','tRNA', 'protein_bind','regulatory', 'misc_recomb','stem_loop','sig_peptide')
        self.tag_qualifiers = ('01-Source','02-Gene','03-mRNA','04-CDS','05-exon','06-intron','07-5UTR','08-CAAT_signal','09-TATA_signal','10-promoter','11-variation','12-polyA_signal','13-polyA_site','14-misc_feature','15-repeat_region','16-precursor_RNA','17-gap','18-prim_transcript','19-mobile_element','20-LTR','21-3UTR','22-misc_difference','23-enhancer','24-mat_peptide','25-rRNA','26-misc_RNA','27-primer_bind','28-tRNA', '29-protein_bind','30-regulatory', '31-misc_recomb','32-stem_loop','33-sig_peptide')
        
        self.listPartial = ['partial CDS', 'partial sequence', 'distal enhancer region', 'promoter region', 'decarboxylase N-terminus', 'flanking region']

                        
    def init_list(self):
        ''' there is sliced - many exons - many translations 
            stores many possible translations '''
        
        self.exon_number_list = []
        
        self.mrna_gene = ''
        self.mrna_product = ''
                
        self.note = ''
        
        self.transl_except = ''
        self.EC_number = ''
        self.locus_tag = ''
        self.old_locus_tag = ''
        self.gene_synonym = ''
        self.standard_name = ''
        self.pseudo = ''
        self.function = ''
        self.inference = ''

                                                
    def outBracket(self, stri):
        return str(stri).replace('[','').replace(']','')
    
    def outBracket_aspas(self, stri):
        stri = str(stri).replace('[','').replace(']','')
        return stri.replace("'",'')

    # sub_features(seq_record, subFeat.featType, dic_sub_feat_exon, feature.sub_features) 
    def proc_sub_features(self, seq_record, featType, feature, showmessage=False):

        '''
            https://github.com/biopython/biopython/blob/master/Bio/SeqFeature.py
            http://biopython.org/DIST/docs/api/Bio.SeqFeature.CompoundLocation-class.html
            
            Rather than sub_features, use a CompoundFeatureLocation
            f.sub_features, f.location should be a CompoundFeatureLocation
              BiopythonDeprecationWarning)
        

        deprecated - 
            if not feature.sub_features:
        '''
    
        if not feature.location:  #.CompoundFeatureLocation:
            return False, '','',[]
    
        ranges = []
        
        if showmessage:
            print '(feature location: ',

            
        if (feature.location.strand == -1): 
            ''' reverse strand '''
            start = feature.location.nofuzzy_end
            end = feature.location.nofuzzy_start
            seq = str(seq_record.seq[end: start:-1])
        else: 
            ''' forward strand '''
            start = feature.location.nofuzzy_start
            end = feature.location.nofuzzy_end-1
            
            ''' ncbi (e.g.) 1:100
                python      0:99
                [start-1: end]
            '''
            seq = str(seq_record.seq[start: end+1])
                    
        ranges.append([start, end])
        
        try:
            gene = feature.qualifiers['gene']
        except:
            gene = ''
        
        try:
            number = feature.qualifiers['number']
        except:
            number = ''
               
        return seq, ranges, gene, number

      
    def proc_sub_features_cds(self, seq_record, feature, showmessage=False):
        if not feature.location:  #.CompoundFeatureLocation:
            return False, '','',[]
    
        seq_list =[]
        ranges = []
            
        for pos in feature.location.parts:
            stri = str(pos.start).replace('<','')
            start = int(stri)
            
            stri = str(pos.end).replace('>','')
            try:
                end = int(stri)
            except:
                try:
                    if end.find('<') >= 0:
                        end = end.replace('<','')
                        end = int(stri) - 1
                    elif end.find('>') >= 0:
                        end = end.replace('>','')
                        end = int(stri) + 1   
                    else:
                        end = 0
                except:
                    end = 0        

            ranges.append([pos.start, pos.end])
            seq_list.append(str(seq_record.seq[start: end+1]))
                
            
        if showmessage:
            print ' - end CompoundFeatureLocation)',
            
        return seq_list, ranges
          

    def do_translation(self, seq, codon_start, showmessage=False):
        seq = seq[codon_start-1:]
        
        dif = len(seq)%3
        if dif>0:
            seq = seq[0: len(seq)-dif]
            
        translated = str(Seq.translate( Seq(seq)))
        
        if showmessage:
            print 'do_translation seq:  %s \n prot: %s \n'%(seq, translated)

        return translated

    def proc_qualifiers(self, featType, feature, int_dic, gi, showmessage=False):
        if not feature.qualifiers:
            return False

        kind = featType + '-quali'
        
        if kind not in int_dic.keys():
            int_dic[kind] = {}
            
        int_dic_qual = int_dic[kind]

        self.loop_features(feature, int_dic_qual)
        
        if showmessage:
            print '  (quali from', featType,
                
        '''
        if featType == 'source':
            self.loop_features_Source(feature, int_dic_qual)
        elif featType == 'CDS':
            self.loop_features_CDS(feature, int_dic_qual)

        elif featType == 'gene':
            self.loop_features_Gene(feature, int_dic_qual)
        elif featType == 'mRNA':
            self.loop_features_mRNA(feature, int_dic_qual)
        elif featType == 'exon':
            self.loop_features_Exon(feature, int_dic_qual)
        elif featType == 'intron':
            self.loop_features_Intron(feature, int_dic_qual)
        else:
            print '\n------------------------------'
            print '>>> featType not found:', featType
            print '------------------------------\n'
            return
        '''

            
        return True
                    
                    
    def loop_features_Source(self, feature, int_dic_qual):
        '''   source: db_xref, ecotype, mol_type, note, organism, isolate '''
        for key in feature.qualifiers.keys():
            ''' gene, note, codon_start, product, protein_id, db_xref, translation, subfeature '''
            int_dic_qual[key] = feature.qualifiers[key]
            
            # print 'source-key', key
            
            if key == 'db_xref':
                self.db_xref = int_dic_qual[key][0]
            elif key == 'ecotype':
                self.ecotype = feature.qualifiers[key][0]
            elif key == 'mol_type':
                self.mol_type = feature.qualifiers[key][0]
            elif key == 'note':
                self.note = feature.qualifiers[key][0]
            elif key == 'organism':
                self.organism_name = feature.qualifiers[key][0]
            elif key == 'isolate':
                self.isolate = feature.qualifiers[key][0]
            elif key == 'clone':
                self.isolate = feature.qualifiers[key][0]
            elif key == 'cultivar':
                self.cultivar = feature.qualifiers[key][0]
            else:
                print '\n=================================='
                print '  <<<<error>>> feature not found %s'%(key)
                print '\==================================\n'


  
    def loop_features(self, feature, int_dic_qual):
        '''   gene: gene '''
        
        for key in feature.qualifiers.keys():
            if key not in int_dic_qual.keys():
                int_dic_qual[key] = []
                
            int_dic_qual[key].append(feature.qualifiers[key][0])
            
                            
    def loop_features_Gene(self, feature, int_dic_qual):
        '''   gene: gene '''
        
        for key in feature.qualifiers.keys():
            int_dic_qual[key] = feature.qualifiers[key][0]
            
            if key == 'gene':
                self.this_gene = int_dic_qual[key]
            else:
                print '\n---------------------------------------'
                print '  <<<<error>>> Gene not found %s'%(key)
                print '---------------------------------------\n'
                return

    def loop_features_mRNA(self, feature, int_dic_qual):
        '''   mRNA: gene, self.product_list '''
        
        for key in feature.qualifiers.keys():
            int_dic_qual[key] = feature.qualifiers[key][0]
            
            if key == 'gene':
                self.mrna_gene = int_dic_qual[key]  
            elif key == 'product':
                self.mrna_product = int_dic_qual[key]   
            else:
                print '\n---------------------------------------'
                print '  <<<<error>>> mRNA not found %s'%(key)
                print '---------------------------------------\n'
                return

         
    def loop_features_Exon(self, feature, int_dic_qual):
        '''   exon:  '''
        for key in feature.qualifiers.keys():
            int_dic_qual[key] = feature.qualifiers[key]
            
            if key == 'gene':
                self.exon_gene = int_dic_qual[key]  
            elif key == 'number':
                #print int_dic_qual[key]
                self.exon_number_list.append(int_dic_qual[key][0])   
            else:
                print '\n---------------------------------------'
                print '  <<<<error>>> Exon not found %s'%(key)
                print '---------------------------------------\n'
                return

    def loop_features_Intron(self, feature, int_dic_qual):
        '''   Intron:  '''
        for key in feature.qualifiers.keys():
            int_dic_qual[key] = feature.qualifiers[key]
            
            if key == 'gene':
                self.exon_gene = int_dic_qual[key]  
            elif key == 'number':
                #print int_dic_qual[key]
                self.exon_number_list.append(int_dic_qual[key][0])   
            else:
                print '\n---------------------------------------'
                print '  <<<<error>>> Intron not found %s'%(key)
                print '---------------------------------------\n'
                return

                                   
    def read_and_parse(self, desk, cntShow=-1, maxiSeq=100000, minSeq = 50, stop=False, showmessage=False):
        cnt, numOk, numError, numPartial = 0,0,0,0
        self.dic_sequence = {}
        
        if not os.path.isfile(self.filenameGBK):
            desk.showmsg_obs('File %s was not found'%(self.filenameGBK))
            return False
        
        for seq_record in SeqIO.parse(open(self.filenameGBK,"r"), "genbank") :
            cnt += 1
            self.dic_sequence[cnt] = {}
            dic_feature_loc = self.dic_sequence[cnt]
            
            msg = '>>>> %i  gi: '%(cnt)
          
            dic_feature_loc['Ok'] = False
            dic_feature_loc['is_partial'] = False
            
            for term in self.listPartial:
                dic_feature_loc[term] = False
                              
            '''
            if cnt == cntShow and cntShow >= 0:
                showmessage = True
            else:
                if cnt > cntShow and cntShow >= 0 and stop:
                    print '---------- EXIT FOR TEST --------------'
                    exit()
                
                showmessage = False
            '''
            
            dic_annotations = {}
            self.init_list()

            accession = ''
            is_partial = False
            gi = ''
            self.species = ''
            self.title = seq_record.description
            
            for term in self.listPartial:
                if self.title.find(term) >= 0:
                    dic_feature_loc[term] = True
                    is_partial = True

            if is_partial:
                numPartial += 1

            for k3 in seq_record.annotations.keys():
                dic_annotations[k3] = self.outBracket_aspas(seq_record.annotations[k3])
                # print k3, dic_annotations[k3]
                
                if k3 == 'accessions':
                    accession = dic_annotations[k3]
                elif k3 == 'gi':
                    gi = dic_annotations[k3]
                    gi = gi.replace(' ','_')
                    gi = gi.replace('-','_')
                    gi = gi.replace(':','')
                    gi = gi.replace('.','')
                    gi = gi.replace(',','')
                    gi = gi.replace(';','')
                    gi = gi.replace('/','')
                    gi = gi.replace('\\','')
                    
                    desk.showmsg_obs(msg + gi)
                elif k3 == 'organism':
                    #if "nebulosa" in dic_annotations[k3]:
                    #    pass
                    stri = self.organism + '_'
                    self.species = dic_annotations[k3].replace(stri,'')
                    self.species = self.species.replace(self.organism,'').strip()
                    self.species = self.species.replace(" ","_")
                    self.species = self.species.replace('-','_')
                    self.species = self.species.replace(':','')
                    self.species = self.species.replace('.','')
                    self.species = self.species.replace(',','')
                    self.species = self.species.replace(';','')
                    self.species = self.species.replace('/','')
                    self.species = self.species.replace('\\','')
                    pass

            if showmessage:
                print '\n=================================' 
                print '%i) gi:%s  accession:%s  organism:%s %s' % (cnt, gi, accession, self.organism, self.species)
             
            dic_feature_loc['accession'] = accession
            dic_feature_loc['gi'] = gi
            dic_feature_loc['organism'] = self.organism
            dic_feature_loc['species'] = self.species
            dic_feature_loc['title'] = self.title
            dic_feature_loc['is_partial'] = is_partial
            dic_feature_loc['note'] = ''
            dic_feature_loc['transl_table'] = None
            dic_feature_loc['transl_except'] = -1
            dic_feature_loc['EC_number'] = ''
            dic_feature_loc['locus_tag'] = ''
            dic_feature_loc['old_locus_tag'] = ''
            dic_feature_loc['gene_synonym'] = ''
            dic_feature_loc['standard_name'] = ''
            dic_feature_loc['pseudo'] = ''
            dic_feature_loc['function'] = ''
            dic_feature_loc['inference'] = ''
            
            to_big_to_little = False

            for feature in seq_record.features:
                featType = self.outBracket_aspas(feature.type)
                try:
                    strand = int(feature.strand)
                except:
                    strand = 1

                if feature.location:
                    loc = feature.location
                    
                    this_ini = loc.nofuzzy_start
                    this_end = loc.nofuzzy_end-1

                    loc_limits = str(this_ini+1) + ':' + str(this_end+1)
                        
                    if showmessage:
                        print '%s(%s)' % (feature.type, strand)
                            
                    if featType == 'source':
                        if len(seq_record.seq) > maxiSeq:
                            print '---------- problems - to BIG %i --------------' % (len(seq_record.seq))
                            print str(seq_record.seq)[:200]
                            print '----------------------------------------------'
                            del(self.dic_sequence[cnt])
                            to_big_to_little = True
                            break
                        elif len(seq_record.seq) < minSeq:
                            print '---------- problems - to LITTLE %i < %i --------------' % (len(seq_record.seq), minSeq)
                            print str(seq_record.seq)
                            print '----------------------------------------------'
                            del(self.dic_sequence[cnt])
                            to_big_to_little = True
                            break
                        
                        if showmessage:
                            print '%s  (%s)' % (feature.type, strand), 'loc_limits:', loc_limits
                        type_num = '01-Source'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]
        
                        #int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        #self.seqCDS = int_dic['seq']
                        int_dic['location'] = loc_limits
                        int_dic['strand'] = strand
                        
                        if showmessage:
                            print '--> source %s ' % (str(int_dic))
                            print dic_feature_loc.keys()

                        self.proc_qualifiers(featType, feature, int_dic, gi, showmessage)
                        ''' source all sequence'''
                        seq, range_lst, gene, number = self.proc_sub_features(seq_record, featType, feature, showmessage=showmessage)

                        int_dic['seq_list'] = [gene, number, seq, range_lst]


                    elif featType == 'gene':
                        type_num = '02-Gene'
                        
                        if type_num not in dic_feature_loc.keys():
                            dic_feature_loc[type_num] = {}
                            int_dic = dic_feature_loc[type_num]
                            
                            int_dic['seq_list'] = []
                            int_dic['list_location'] = []
                            
                        else:
                            int_dic = dic_feature_loc[type_num]
                        

                        self.proc_qualifiers(featType, feature, int_dic, gi, showmessage=False)
                                
                        ''' gene - all gene '''
                        seq, range_lst, gene, number = self.proc_sub_features(seq_record, featType, feature, showmessage=showmessage)

                        int_dic['seq_list'].append([gene, number, seq, range_lst])

                            
                        if showmessage:
                            print '--> gene %s ' % (str(dic_feature_loc[type_num]))
                            print dic_feature_loc.keys()

                    elif featType == 'mRNA':
                        type_num = '03-mRNA'
                        
                        if type_num not in dic_feature_loc.keys():
                            dic_feature_loc[type_num] = {}
                            int_dic = dic_feature_loc[type_num]

                            int_dic['seq_list']  = []
                            int_dic['list_location']  = []                         
                        else:
                            int_dic = dic_feature_loc[type_num]

                        self.proc_qualifiers(featType, feature, int_dic, gi, showmessage)

                        ''' mRNA '''
                        seq, range_lst, gene, number = self.proc_sub_features(seq_record, featType, feature, showmessage=showmessage)
                        int_dic['seq_list'].append([gene, number, seq, range_lst])

                        if showmessage:
                            print '--> mRNA %s ' % (str(dic_feature_loc[type_num]))
                            print dic_feature_loc.keys()

                    elif featType == 'CDS':

                        type_num = '04-CDS'
                        if type_num not in dic_feature_loc.keys():                        
                            dic_feature_loc[type_num] = {}
                            int_dic = dic_feature_loc[type_num]
                            
                            int_dic['seq_list'] = []
                            int_dic['location_list'] = []           
                        else:
                            int_dic = dic_feature_loc[type_num]

                        self.proc_qualifiers(featType, feature, int_dic, gi, showmessage)

                        ''' CDS '''
                        seq_list, range_lst = self.proc_sub_features_cds(seq_record, feature, showmessage=showmessage)
                        
                        int_dic['seq_list'].append(seq_list)
                        int_dic['location_list'].append(range_lst)
                                                    

                        # print '\n   ........................................'
                        # description = organism + '_' + self.product_list
                        # seq = SeqRecord(self.seqExon_CDS, id=gi, name=organism, description=description)
                        if showmessage:
                            print '---> CDS %s %s ' % (type_num, str(dic_feature_loc[type_num]))
                            print dic_feature_loc.keys()
                        
                    elif featType == 'exon':
                        type_num = '05-exon'
                        
                        if type_num not in dic_feature_loc.keys():
                            dic_feature_loc[type_num] = {}
                            int_dic = dic_feature_loc[type_num]
                            
                            int_dic['seq_list'] = []
                            int_dic['list_location'] = []
                            
                        else:
                            int_dic = dic_feature_loc[type_num]                            

                        self.proc_qualifiers(featType, feature, int_dic, gi, showmessage)

        
                        ''' part of exons '''
                        seq, range_lst, gene, number = self.proc_sub_features(seq_record, featType, feature, showmessage=showmessage)

                        int_dic['seq_list'].append([gene, number, seq, range_lst])


                    elif featType == 'intron':
                        type_num = '06-intron'
                        
                        if type_num not in dic_feature_loc.keys():
                            dic_feature_loc[type_num] = {}
                            int_dic = dic_feature_loc[type_num]
                            
                            int_dic['seq_list'] = []
                            int_dic['list_location'] = []                         
                        else:
                            int_dic = dic_feature_loc[type_num]


                        self.proc_qualifiers(featType, feature, int_dic, gi, showmessage)

                        seq, range_lst, gene, number = self.proc_sub_features(seq_record, featType, feature, showmessage=showmessage)
                        
                        int_dic['seq_list'].append([gene, number, seq, range_lst])
                        
        
                    elif featType == '5UTR':
                        type_num = '07-5UTR'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits

        
                    elif featType == 'CAAT_signal':
                        type_num = '08-CAAT_signal'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits
        
                    elif featType == 'TATA_signal':
                        type_num = '09-TATA_signal'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits
        
                    elif featType == 'promoter':
                        type_num = '10-promoter'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits
        
                    elif featType == 'variation':
                        type_num = '11-variation'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits
        
                    elif featType == 'polyA_signal':
                        type_num = '12-polyA_signal'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits

                    elif featType == 'polyA_site':
                        type_num = '13-polyA_site'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits
                        
                    elif featType == 'misc_feature':
                        type_num = '14-misc_feature'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits
                        
                    elif featType == 'repeat_region':
                        type_num = '15-repeat_region'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits

                    elif featType == 'precursor_RNA':
                        type_num = '16-precursor_RNA'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits
                    elif featType == 'gap':
                        type_num = '17-gap'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits
                    elif featType == 'prim_transcript':
                        type_num = '18-prim_transcript'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits
                        
                    elif featType == 'mobile_element':
                        type_num = '19-mobile_element'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits
                        
                    elif featType == 'LTR':
                        type_num = '20-LTR'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits
                        
                    elif featType == '3UTR':
                        type_num = '21-3UTR'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits
                        
                    elif featType == 'misc_difference':
                        type_num = '22-misc_difference'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits

                    elif featType == 'enhancer':
                        type_num = '23-enhancer'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits

                    elif featType == 'mat_peptide':
                        type_num = '24-mat_peptide'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits
                    
                    elif featType == 'rRNA':
                        type_num = '25-rRNA'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits

                    elif featType == 'misc_RNA':
                        type_num = '26-misc_RNA'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits

                    elif featType == 'primer_bind':
                        type_num = '27-primer_bind'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits

                    elif featType == 'tRNA':
                        type_num = '28-tRNA'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits

                    elif featType == 'protein_bind':
                        type_num = '29-protein_bind'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits


                    elif featType == 'regulatory':
                        type_num = '30-regulatory'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits

                    elif featType == 'misc_recomb':
                        type_num = '31-misc_recomb'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits


                    elif featType == 'stem_loop':
                        type_num = '32-stem_loop'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits

                    elif featType == 'sig_peptide':
                        type_num = '33-sig_peptide'
                        dic_feature_loc[type_num] = {}
                        int_dic = dic_feature_loc[type_num]

                        int_dic['seq'] = str(seq_record.seq[this_ini: this_end+1])
                        int_dic['location'] = loc_limits

                    elif featType == 'old_sequence':
                        pass
                    
                    elif featType == 'STS':
                        pass

                    
                    else:
                        print '\n---------------------------------------'
                        print ">>>> new key int_dic ??? for '%s'"%(featType)
                        print '---------------------------------------\n'



            #-------------------------------------------------------------------------
            
            num = 0

            try:
                codon_start_list = dic_feature_loc['04-CDS']['CDS-quali']['codon_start']
            except:
                codon_start_list = []

            try:
                gene_list = dic_feature_loc['04-CDS']['CDS-quali']['gene']
            except:
                gene_list = []

            try:
                product_list = dic_feature_loc['04-CDS']['CDS-quali']['product']
            except:
                product_list = []
                
            try:
                protein_id_list = dic_feature_loc['04-CDS']['CDS-quali']['protein_id']
            except:
                protein_id_list = []

                
            if gene_list:
                gene = gene_list[num]
            else:
                gene = self.gene
                
            if product_list:
                product = product_list[num]
            else:
                product = 'unknown'
                
                
            if protein_id_list:
                protein_id = protein_id_list[num]
            else:
                protein_id = 'unknown'   

            dic_feature_loc['gene'] = gene
            dic_feature_loc['product'] = product
            dic_feature_loc['protein_id'] = protein_id
            
            '''===== analysis translation ======='''
                                
            try:
                codon_start = int(codon_start_list[num])
            except:
                codon_start = 1
                
                            
            try:
                listTranslations = dic_feature_loc['04-CDS']['CDS-quali']['translation']
                seqTranslation = listTranslations[num]
            except:
                seqTranslation = ''

                
            try:
                # listListCDS = dic_feature_loc['04-CDS']['seq_list']
                ''' gene, number, seq, range_lst '''
                listListExon = dic_feature_loc['05-exon']['seq_list']
                seqExon = ''
                for mat in listListExon:
                    seqExon += mat[2]
            except:
                try:
                    listCDS = dic_feature_loc['04-CDS']['seq_list'][num]
                    
                    seqExon = ''
                    for seq in listCDS:
                        seqExon += seq
                except:
                    seqExon = ''

                    
            if seqExon <> '':    
                translated = self.do_translation(seqExon, codon_start, showmessage=showmessage)
            else:
                translated = ''
                
            ''' there may be many .... getting the first '''

              
            if showmessage:    
                print 'do_translation'

                print ' CDS  --> ', seqTranslation
                print ' EXON --> ', translated

            if to_big_to_little:
                print 'To big or to little'
                numError+=1
                continue
                
            if seqTranslation and translated:
                lenTransl = len(translated)
                lenAA = len(seqTranslation)
                
                ''' get the minor  = both must be equal in length'''
                if lenAA <> lenTransl:
                    if lenAA < lenTransl:
                        print '>>> correction !!! lenAA %i < lenTransl %i'%(lenAA, lenTransl)
                        lenTransl = lenAA
                        seqTranslation = seqTranslation[:lenAA]
                    else:
                        print '>>> correction !!! lenAA %i > lenTransl %i'%(lenAA, lenTransl)
                        lenAA = lenTransl
                        translated = translated[:lenTransl]

                    
                    # print ' --> ', seqTranslation 
            else:
                lenAA = 0


            if  seqTranslation  == '':
                if not is_partial:
                    numError+=1
                    continue
                    
            elif seqTranslation == translated:
                numOk += 1
                dic_feature_loc['Ok'] = True
                
                if showmessage:
                    print '------ All right -----'
            else:
                ''' if errors '''
                if showmessage:
                    print '\n Descri:', dic_feature_loc['title']
                    print ' Source:', dic_feature_loc['01-Source']['seq']
                    print '   range', range_lst
                    print ' This AA', seqTranslation
                    # print ' Tr.NCBI', cds_translationList[counter_cds]
                    
                    print '\ncds_translation         %s' %(seqTranslation)
                    print 'seqTranslation 1st %s' %(seqTranslation)
                    

                continue
                                  
        # print '\n-------------------------------------------------------------'              
        # print "%i seqs, OK %i, Error %i, Partial %i, Error=%5.2f%%  Empty=%5.2f%%  "%(cnt, numOk, numError, numPartial, 100.*numError/cnt, 100.*((cnt-numOk)/float(cnt)))
        # print '-------------------------------------------------------------\n\n'              
        return True
    
    def show_data(self, cntShow=-1, stop=False,showmessage=False):
        numError, withoutTranslation, numPartial, numOk = 0, 0, 0, 0
        
        multi_CDS = False
        multi=[]
        
        for cnt in self.dic_sequence.keys():
            '''
            if cnt == cntShow and cntShow >= 0:
                showmessage = True
            else:
                if cnt > cntShow and cntShow >= 0 and stop:
                    print '---------- EXIT FOR TEST --------------'
                    exit()
                showmessage = False
            '''
            
            dic_feature_loc = self.dic_sequence[cnt]
            
            num = 0
            try:
                seqTranslation = dic_feature_loc['04-CDS']['CDS-quali']['translation'][num]
            except:
                seqTranslation = ''
            
            try:
                numCDS = len(dic_feature_loc['04-CDS']['seq_list'])
                if numCDS > 1:
                    multi_CDS = True
                    multi.append([dic_feature_loc['gi'], dic_feature_loc['accession'], numCDS])
                else:
                    multi_CDS = False
            except:
                multi_CDS = False
                
                
            hasCDS = (dic_feature_loc['Ok'] == True)
            if hasCDS:
                numOk += 1
                
                if showmessage:
                    print cnt, ') OK %s' % (dic_feature_loc['gi'])
                    
                continue
            else:
                if dic_feature_loc['is_partial'] == True:
                    numPartial += 1
                    continue
                
                if seqTranslation == '':
                    withoutTranslation += 1
                    print '%i) ??? %s has no translation'%(cnt, dic_feature_loc['gi'])
            
            for key in dic_feature_loc:
                if key == 'gi' or key == 'accessions':
                    continue
                
                if type(dic_feature_loc[key]) <> type({}):
                    if showmessage:
                        print '   .', key, dic_feature_loc[key]  # rever
                    continue

                dic = dic_feature_loc[key]
                if showmessage:
                    print '   .', key 
                for key2 in dic.keys():
                    if not dic[key2]:
                        continue
                    
                    if showmessage:
                        print '      .', key2,  # type(dic[key2])
                    
                    if type(dic[key2]) == type('a'):
                        if showmessage:
                            print ' - ', dic[key2]
                        '''
                        if len(dic[key2]) > 25:
                            print ' - ', dic[key2][:25], '...'
                        else:
                            print ' - ', dic[key2]
                        '''
                    elif type(dic[key2]) == type(list()):
                        if showmessage:
                            print ' - ', dic[key2]
                    elif type(dic[key2]) == type(True):
                        if showmessage:
                            print ' - ', dic[key2]
                    elif type(dic[key2]) == type(2):
                        if showmessage:
                            print ' - ', dic[key2]
                    else:
                        if showmessage:
                            for key3 in dic[key2]:
                                print '         .', key3, dic[key2][key3]
                             
            if numError == 3:
                print '\n---------------------------------------'
                print '>>>>> error:  STOP WITH 3 ERRORS '
                print '---------------------------------------\n'
                return
            
            '''---- do you want control translations ??? ---
            if withoutTranslation == 300 :
                print '\n---------------------------------------'
                print '>>>>> error: STOP WITH 300 NO Translations '
                print '---------------------------------------\n'
                return
            -----------------------------------------------'''

        print '\n-------------------------------------------------------------'              

        if cnt > 0:             
            print "%i seqs, OK %i, Error %i, Partial %i, NoTransl %i, Error=%5.2f%%  Problems=%5.2f%%  "%(cnt, numOk, numError, numPartial, withoutTranslation, 100.*numError/cnt, 100.*((cnt-numOk)/float(cnt)))
        else:
            print "0 seqs"
        print '-------------------------------------------------------------'              

        if multi_CDS:
            print 'Multi CDS'
            for elem in multi:
                print 'gi: %s   accession: %s  num multiple CDS: %i' % (elem[0], elem[1], elem[2])
        else:
            print 'There is no Multi CDS'
        print '-------------------------------------------------------------'              

    def save_data_cds_exon_intron_prot(self, desk, includePartial, showmessage=False):
        dic_species_params = {}
        dic_record_Gene = {}
        dic_record_CDS = {}
        dic_record_Exon = {}
        dic_record_Intron = {}
        dic_record_Protein = {}

        gi_list = []
        count_saves = 0
        count_already_exists = 0
        
        for cnt in self.dic_sequence.keys():
            
            dic_feature_loc = self.dic_sequence[cnt]        
            
            gi = dic_feature_loc['gi']
            
            if gi in gi_list:
                count_already_exists += 1
                print 'This %s already exists.'%(gi)
                continue
            
            gi_list.append(gi)
            '''
            if cnt == 81:
                print cnt, gi, dic_feature_loc.keys()
            '''
            try:
                species = dic_feature_loc['species']
            except:
                species = 'unknown'            
 
            try:
                gene = dic_feature_loc['gene']
            except:
                gene = 'unknown' 
                
            try:
                product = dic_feature_loc['product']
            except:
                product = 'unknown' 
                
            try:
                _ = dic_feature_loc['protein_id']
            except:
                _ = 'unknown' 
                
            
            if species not in dic_species_params.keys():
                dic_species_params[species] = ['x',0,0,0,0,0,0]
                
            
            desc = dic_feature_loc['organism'] + '|' + dic_feature_loc['species'] + '|' + gene
            desc = desc.replace(' ', '-');
            
            desc += '|Prod|' + product # '||' + dic_feature_loc['title']
            desc = desc.replace(' ', '_')
            desc = desc.replace('-', '_')
            desc = desc.replace('\\', '')
            desc = desc.replace('/', '')
            desc = desc.replace(':', '')
            desc = desc.replace('.', '')
            desc = desc.replace(',', '')
            desc = desc.replace(';', '')

            gi = gi + '_' + dic_feature_loc['accession'] + '_' + desc
            desc = ''
        
            if gi.lower().find('synthetic') >= 0:
                print '  !!! Do not accept synthetic: %s'%(gi)
                continue
            
            
            try:
                is_partial = dic_feature_loc['is_partial']
            except:
                is_partial = True

            ''' getting the 1st '''
            mySeqCDS = ''
            num = 0
            
            try:
                listCDS = dic_feature_loc['04-CDS']['seq_list']
                
                ''' getting only the first '''
                for seq in listCDS[num]:
                    mySeqCDS += seq                
            except:
                mySeqCDS = ''
                
            try:
                mySeqProt = dic_feature_loc['04-CDS']['CDS-quali']['translation'][num]
            except:
                mySeqProt = ''
                
                
            mySeqExon = ''
            mySeqIntron = ''
            mySeqGene = ''
            try:
                for mat in dic_feature_loc['05-exon']['seq_list']:
                    mySeqExon += mat[2]
            except:
                pass
                
            
            try:   
                for mat in dic_feature_loc['06-intron']['list_seq_intron']:
                    mySeqIntron += mat[2]
            except:
                pass
                         
            try:
                for mat in dic_feature_loc['02-Gene']['seq_list']:
                    mySeqGene += mat[2]
            except:
                pass
                               
            
            '''
            try:
                mySeqMRNA = dic_feature_loc['03-mRNA']['seq'][0]
            except:                                
                mySeqMRNA = ''
            '''

           
            if len(mySeqGene) >= self.cutoffLength and (is_partial == False or includePartial):
                count_saves += 1
                '''
                if gi == '85665937_DQ343316_Capsella-bursa-pastoris_Adh|All|unkAllele|Prod|alcohol-dehydrogenase':
                    print gi
                    print mySeqGene
                    pass
                '''
                if mySeqGene <> '':
                    if species not in dic_record_Gene.keys():
                        dic_record_Gene[species] = []
                    dic_record_Gene[species].append(SeqRecord(Seq(mySeqGene), id=gi, name='', description=''))
                      
                if mySeqCDS <> '':
                    if showmessage:
                        print '%i) cds gi %s, partial %s, exon: %s' % (cnt, gi, str(is_partial), mySeqCDS[1:25]),
                        
                    if species not in dic_record_CDS.keys():
                        dic_record_CDS[species] = []
                    dic_record_CDS[species].append(SeqRecord(Seq(mySeqCDS), id=gi, name='', description=''))
            
                if mySeqExon <> '':
                    if showmessage:
                        print '%i) exon gi %s, partial %s, exon: %s' % (cnt, gi, str(is_partial), mySeqExon[1:25]),

                    if species not in dic_record_Exon.keys():
                        dic_record_Exon[species] = []
                    dic_record_Exon[species].append(SeqRecord(Seq(mySeqExon), id=gi, name='', description=''))
    
                    # seqProt = str(oneSeq.translate(table="Standard", stop_symbol="*", to_stop=False, cds=False))
                    
            
                if mySeqIntron <> '':
                    if species not in dic_record_Intron.keys():
                        dic_record_Intron[species] = []
                    dic_record_Intron[species].append(SeqRecord(Seq(mySeqIntron), id=gi, name="", description=""))
            
                if mySeqProt <> '':
                    if showmessage:
                        print 'prot   AA %s ...' % (mySeqProt[1:25])
                    if species not in dic_record_Protein.keys():
                        dic_record_Protein[species] = []
                    dic_record_Protein[species].append(SeqRecord(Seq(mySeqProt), id=gi, name="", description=""))
        
                     

        mySeq = MySequence(desk, showmessage=False)

        seq_record_Gene = []
        seq_record_CDS = []
        seq_record_Exon = []
        seq_record_Intron = []
        seq_record_Protein = []
        
        
        for species in dic_species_params.keys():

            filename = '%s_%s_%s_%s_%iL' % \
            (self.organism, species, '%s', self.gene_title, self.cutoffLength)
            
            filename = filename.replace(' ','_')
            filename = filename.replace('-','_')
            filename = filename.replace(':','')
            filename = filename.replace('.','')
            filename = filename.replace(',','')
            filename = filename.replace(';','')
            filename = filename.replace('/','')
            filename = filename.replace('\\','')

            if species in dic_record_Gene.keys():
                
                filename_GeneFasta = self.rootFasta + filename%('Gene') + '.fasta'
                
                ''' species and how many seqs there is '''
                mat = dic_species_params[species] 
                mat[1] = len(dic_record_Gene[species])
                mat[2] = len(dic_record_Gene[species])
                mat[3] = len(dic_record_Gene[species])

                mat[4] = len(dic_record_Gene[species][0])
                mat[5] = len(dic_record_Gene[species][0])
                mat[6] = len(dic_record_Gene[species][0])
                
                seq_record_Gene += dic_record_Gene[species]
                mySeq.save_fasta_new(filename_GeneFasta, dic_record_Gene[species])
                
            if species in dic_record_CDS.keys():
                filename_CDSFasta = self.rootFasta + filename%('CDS') + '.fasta'
                seq_record_CDS += dic_record_CDS[species]
                mySeq.save_fasta_new(filename_CDSFasta, dic_record_CDS[species])
                
            if species in dic_record_Exon.keys():
                filename_ExonFasta = self.rootFasta + filename%('Exon') + '.fasta'
                seq_record_Exon += dic_record_Exon[species]
                mySeq.save_fasta_new(filename_ExonFasta, dic_record_Exon[species])
            
            if species in dic_record_Intron.keys():
                filename_IntronFasta = self.rootFasta + filename%('Intron') + '.fasta'
                seq_record_Intron += dic_record_Intron[species]
                mySeq.save_fasta_new(filename_IntronFasta, dic_record_Intron[species])
                
            if species in dic_record_Protein.keys():
                filename_ProtFasta = self.rootFasta + filename%('Prot') + '.fasta'
                seq_record_Protein += dic_record_Protein[species]
                mySeq.save_fasta_new(filename_ProtFasta, dic_record_Protein[species])
        

        ''' saving all sequens in _all_ files '''
        if seq_record_Gene:
            mySeq.save_fasta_new(self.filename_GeneFasta, seq_record_Gene)
        if seq_record_CDS:
            mySeq.save_fasta_new(self.filename_CDSFasta, seq_record_CDS)
        if seq_record_Exon:
            mySeq.save_fasta_new(self.filename_ExonFasta, seq_record_Exon)
        
        
        if seq_record_Intron:
            mySeq.save_fasta_new(self.filename_IntronFasta, seq_record_Intron)
            
        if seq_record_Protein:
            mySeq.save_fasta_new(self.filename_ProtFasta, seq_record_Protein)

        print '-------------------------------------------------------------' 
        desk.showmsg_obs('Saved %i records with at least len(Gene) >= %i and not Partial CDS, repeated %i,  Seqs: Gene  %i, CDS %i, Exon %i, Intron %i, Protein %i'%\
             (count_saves, self.cutoffLength, count_already_exists, len(seq_record_Gene), len(seq_record_CDS), len(seq_record_Exon), len(seq_record_Intron), len(seq_record_Protein)) )
        
        print '-------------------------------------------------------------' 
        
        return dic_species_params


            