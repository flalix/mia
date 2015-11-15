#-*- coding: utf-8 -*-
'''
Created on 07/07/2014
Updated on 21/07/2014
Updated on 27/08/2014
Updated on 01/09/2014
Updated on 22/09/2014
Updated on 24/11/2014
Updated 0n 23/12/2014
Updated on 22/01/2015 - Linux tests
Updated on 24/03/2015 - no log check searchBox
Updated on 07/08/2015 - no more dollar var
Updated on 27/08/2015 - mrbayes and classifier
Updated on 06/10/2015 - get_params()

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica
'''
import os, sys, copy, gc, platform
import ttk, tkMessageBox, tkFont, Tkinter as tk
from Tkinter import END #DISABLED, ACTIVE
from subprocess import call
from pip._vendor.distlib.compat import which
''' # BOTTOM, LEFT, RIGHT, TOP, CENTER, N, NE, E, SE, S, SW, tk.NW, CENTER '''
# from Tkinter import tk.Frame,Radiobutton, tk.Checkbutton, Entry, tk.StringVar, tk.IntVar, tk.BooleanVar
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from classes import mrBayesClass as MB
import classes.Drosophila as dro

import pipe01_Find_ManySpecies_OneGene_GBK as pipe1
import pipe02_Read_GBK_to_Fasta as pipe2
import pipe03_alignment as pipe3
import pipe04_purge_repeated_sequences as pipe4
import pipe05_Change_Ambiguous_Nuc as pipe5
import pipe06_Vertical_MI_HeatMap as pipe6
import pipe07_Horizontal_MI_Distribution as pipe7
import pipe08_JSD as pipe8
import pipe09_Hierarchical_Cluster as pipe9
import pipe10_random_Entropy_Experiment_Shannon_WriteFile as pipe10
import pipe12_mr_bayes_params as pipe12
import pipe13_hmi_classifier as pipe13

tk_root = None
mainFrameWidth = 100

inFrameWidth = mainFrameWidth
lenInputFrame10 = inFrameWidth
lenInputFrame70 = 70

widthSpeciesFrame = 40
widthListSpecies = widthSpeciesFrame 
secondFieldWidth = 5


class Desktop(tk.Frame):
    ''' An example application for TkInter.  Instantiate
        and call the run method to run. 
        http://www.pythoncentral.io/introduction-python-gui-development/
        https://docs.python.org/2/library/tk.html
    '''
    def __init__(self):
        ''' checkbox,  num of seq (ori, min, max)  len (ori,min,max)  species '''
        self.mask = "[%s] %3i/%3i/%3i | %5i/%5i/%5i | %s"
        self.cour10 = tkFont.Font(family='courier', size=10, weight='normal')  # you don't have to use Helvetica or bold, this is just an example
        
        ''' Initialize window using the parent's constructor
            width=mainFrameWidth, height=mainFrameHeight, '''
        tk.Frame.__init__(self, tk_root, relief=tk.SUNKEN,bd=2)


        self.colors = ['aliceblue', 'aqua', 'blue', 'chocolate', 'cyan','darkgreen','green',
                       'gray','red','ivory','indigo','lavender',
                       'lightseagreen', 'olive', 'orange', 'purple','salmon']
       
        self.first_time = True
        self.new_search = False

        '''
        print platform.system()
        '''
        
        self.platform = platform.system() # 'Linux'  Win32
        self.isWindows = (self.platform == 'Windows') 
        
        if self.isWindows:
            print("Is Windows")
        else:
            print("Is Linux ... (Mac not tested)")
        
        # Set the title
        self.tk_root = tk_root
        self.tk_root.columnconfigure(0, weight=1)
        self.tk_root.rowconfigure(0, weight=1)
    
        tk_root.title('Mutual Information Analyzer')
        tk_root.resizable(False,False)
        
        self.speciesListWidth  = 70
        self.speciesListHeight = 30

        self.set_defaults() 
        self.init_vars()
        self.read_default_ini()
        #self.initialize_mia()

        self.nb = ttk.Notebook(tk_root)  # width=inFrameWidth
        
        self.text_list = [ ['NCBI', 'Search Species and Genes in NCBI'], \
                           ['Gbk-Fasta', 'Transform GBK to Fasta'] ,  \
                           ['Alignment', 'Align Muscle or Clustalw'], \
                           ['Purging', 'Purging'],  \
                           ['Consensus', 'Substitute Ambiguous Nucleotides'],\
                           ['VMI', 'VMI - Vertical Mutual Information'], \
                           ['HMI', 'HMI - Horizontal Mutual Information'], \
                           ['JSD', 'Jensen-Shannon divergence'], 
                           ['Cluster', 'Hierarchical Cluster'], 
                           ['Entropy', 'Shannon Entropy'], 
                           ['MrBayes', 'MrBayes params analyzer'], 
                           ['Classifier', 'HMI Classifier'], 
                           ['Root', 'Root']]
        
        numTabs = len(self.text_list)
        
        self.listImage = ['png','tif','jpg']
        self.listHeatMapColors = ['afmhot','cool','Set1','hot','spring','summer','winter','bwr','seismic','Paired','Blues','Reds']


        self.frames = []
        for i in range(numTabs):
            self.frames.append(tk.Frame())  # height=mainFrameHeight, width=inFrameWidth))

        for i in range(numTabs):
            self.nb.add(self.frames[i], text=self.text_list[i][0])
            self.start(self.frames[i],i)

        
        self.nb.grid(row=0, column=0, sticky='NWE')

        self.tab_caption = 'NCBI'

        self.drawButtons(row=2)
        ''' must be the last, to calc the height '''
        self.drawText(row=1)
        
        self.drawSpeciesList()

        self.build_listSpecies()
        self.filter_listSpecies(mode=self.rb_seq_filter_var.get(), cutoff=self.cutoffNumSeq_var.get(), showmessage=True)

        self.nb.bind_all('<Key>', self.on_change_fields)        
        self.nb.bind_all("<<NotebookTabChanged>>", self.tabChangedEvent)


        if self.organism == "Drosophila":
            self.dr = dro.Drosophila()
        else:
            self.dr = None

        self.mb = MB.mrBayesClass(self)
        
        # attach cursor to this widget
        self.cursor = tk_root.cget("cursor")
         
    def warning(self, stri):
        tkMessageBox.showinfo('warning', stri)
        
        
    def showmsg_obs(self, msg_obs, same_line = False):
        if type(msg_obs) != type("a"):
            msg_obs = str(msg_obs)
            
        if same_line:
            print msg_obs,
        else:
            print msg_obs
        
        self.msg_obs = msg_obs
        self.msg_obs_var.set(msg_obs)
        tk_root.update()


    def init_vars(self):
        #self.labelroot_var = tk.StringVar()
        self.dna_prot_var = tk.StringVar()
        self.organism_var = tk.StringVar()
        self.gene_var = tk.StringVar() 
        self.species_var = tk.StringVar()
        self.cutoffLength_var = tk.IntVar()
        self.cutoffNumSeq_var = tk.IntVar()
        self.numOfLetters_var = tk.IntVar()
        self.rootFasta_var = tk.StringVar()
        self.rootImage_var = tk.StringVar()
        self.rootTable_var = tk.StringVar()
        self.rootTree_var = tk.StringVar()
        self.rootEntropy_var = tk.StringVar()
        self.seqType_var = tk.StringVar()
        self.showGraph_var = tk.BooleanVar()
        self.saveGraph_var = tk.BooleanVar()
        self.de_novo_var = tk.BooleanVar()

        self.geneList_var = tk.StringVar()
        self.specieList_var = tk.StringVar()
        self.search_org_gene_var = tk.StringVar()
        self.db_var = tk.StringVar()
        
        self.partialCDS_var = tk.StringVar()
        self.cds_var = tk.StringVar()
        self.mrna_var = tk.StringVar()
        self.completeSeq_var = tk.StringVar()
        self.completeGen_var = tk.StringVar()
        self.chromosome_var = tk.StringVar()
        self.shotGun_var = tk.StringVar()
        self.contig_var = tk.StringVar()
        self.superContig_var = tk.StringVar()
        
        self.logGbkFilename_var = tk.StringVar()
        self.retmax_var = tk.IntVar()
        self.email_var = tk.StringVar()
         
        self.maxiSeq_var = tk.IntVar()
        self.minSeq_var  = tk.IntVar()
        
        self.stop_var = tk.BooleanVar()
        self.showmessage_var = tk.BooleanVar()
        self.confirm_gene_var = tk.BooleanVar()

        self.listBadTerms_var = tk.StringVar()
        
        self.output_filename_gbk_var = tk.StringVar()
        
        self.input_filename_align_var = tk.StringVar()
        self.output_filename_purged_var = tk.StringVar()
        self.output_filename_fasta_var = tk.StringVar()
        self.aligned_filename_fasta_var = tk.StringVar()
        self.output_filename_consensus_var = tk.StringVar()
        
        self.output_hmi_var = tk.StringVar()
        self.output_vmi_var = tk.StringVar()
        
        self.input_jsd_var = tk.StringVar()

        self.numOfLetters_var = tk.IntVar()
        self.frame_var = tk.IntVar()
        self.isLog_var = tk.BooleanVar()
        self.offset_var = tk.IntVar()

        self.msg_obs_var = tk.StringVar()
        self.totList_var = tk.IntVar()
        
        self.totList_var = tk.IntVar()
        self.totSeqs_var = tk.IntVar()
        
        self.rb_seq_filter_var = tk.StringVar()
        self.minmax_var = tk.StringVar()
        self.each_all_var = tk.StringVar()
        
        self.rb_inv_filter_var  = tk.StringVar()
                
        self.vert_horiz_var = tk.StringVar()
        
        self.saveData_var = tk.BooleanVar()

        self.numOfExperiments_var = tk.IntVar()
        self.sartAt_var = tk.IntVar()
        self.lengthSim_var = tk.IntVar()
        
        self.shannon_filename_var = tk.StringVar()
        
        self.title_var = tk.StringVar()
                
        self.cluster_input_filename_var = tk.StringVar()
        self.cluster_method_var = tk.StringVar()

        self.withCorrection_var = tk.BooleanVar()
        self.mnat_var = tk.BooleanVar()
        self.numOfSEs_var = tk.IntVar()
        self.norm_var = tk.BooleanVar()
        
        self.imgType_var = tk.StringVar()
        self.dpi_var = tk.IntVar()
        self.colorThreshold_var = tk.DoubleVar()
        
        self.heatmap_color_var = tk.StringVar()
        self.heatmap_ceil_value_var = tk.DoubleVar()
        self.heatmap_ceil_var = tk.BooleanVar()
        self.is3D_var = tk.BooleanVar()
       
                
        self.alignment_var = tk.StringVar()
        self.minVertCutoff_var = tk.IntVar()
        self.maxVertCutoff_var = tk.IntVar()
        self.horizCutoff_var = tk.IntVar()

        self.realing_each_species_var = tk.BooleanVar()
        self.realing_all_var = tk.BooleanVar()

        self.start_stop_codon_var = tk.BooleanVar()
        
        ''' mrBayes '''
        self.mrBayes_path_var = tk.StringVar()
   
        self.piA_var = tk.BooleanVar()
        self.piC_var = tk.BooleanVar()
        self.piG_var = tk.BooleanVar()
        self.piT_var = tk.BooleanVar()
        self.rAC_var = tk.BooleanVar()
        self.rAG_var = tk.BooleanVar()
        self.rAT_var = tk.BooleanVar()
        self.rCG_var = tk.BooleanVar()
        self.rCT_var = tk.BooleanVar()
        self.rGT_var = tk.BooleanVar()
        self.LnL_var = tk.BooleanVar()
        self.LnPr_var = tk.BooleanVar()
        self.TL_var = tk.BooleanVar()
        self.alpha_var = tk.BooleanVar()
        self.off_on_var = tk.BooleanVar()
        self.on_off_var = tk.BooleanVar()
        self.pinvar_var = tk.BooleanVar()
        self.LnL_HMean_var = tk.BooleanVar()
        
        ''' classifier '''
        self.NminSamples_var = tk.IntVar()
        self.numSimLoops_var = tk.IntVar()
        
        self.figtree_path_var  = tk.StringVar() 
        self.mrBayes_replace01_var = tk.StringVar() 
        self.mrBayes_replace02_var = tk.StringVar() 
        self.mrBayes_replace03_var = tk.StringVar() 
        self.mrBayes_replace04_var = tk.StringVar() 
        self.string_kill_var = tk.StringVar()     # '_Adh_Prod_
        self.colaps_tree_var = tk.BooleanVar()
        self.aligned_consensus_var = tk.StringVar() 
        self.standard_covarion_var = tk.StringVar() 
        
        self.rand_method_var = tk.StringVar() 
        
        self.mrBayes_fn01_var = tk.StringVar()         
        self.mrBayes_fn02_var = tk.StringVar()         
        self.mrBayes_fn03_var = tk.StringVar()         
        self.mrBayes_fn04_var = tk.StringVar()

        self.mrBayes_fn01_cov_var = tk.StringVar()         
        self.mrBayes_fn02_cov_var = tk.StringVar()         
        self.mrBayes_fn03_cov_var = tk.StringVar()         
        self.mrBayes_fn04_cov_var = tk.StringVar()
      
        self.burnin_var = tk.DoubleVar()      

        
    def start(self, tab, num):
        self.tab_caption = self.text_list[num][0]
        
        if self.tab_caption == 'NCBI':
            self.form_NCBI(num, tab)
        elif self.tab_caption == 'Gbk-Fasta':
            self.form_GBK_Fasta(num, tab)
        elif self.tab_caption == 'Alignment':
            self.form_aligment(num, tab)
        elif self.tab_caption == 'Purging':
            self.form_Purging(num, tab)
        elif self.tab_caption == 'Consensus':
            self.form_Consensus(num, tab)
        elif self.tab_caption == 'HMI':
            self.form_HMI(num, tab)
        elif self.tab_caption == 'VMI':
            self.form_VMI(num, tab)
        elif self.tab_caption == 'JSD':
            self.form_JSD(num, tab)
        elif self.tab_caption == 'Cluster':
            self.form_cluster(num, tab)
        elif self.tab_caption == 'Entropy':
            self.form_entropy(num, tab)
        elif self.tab_caption == 'Root':    
            self.form_root(num, tab)
        elif self.tab_caption == 'MrBayes':    
            self.form_mr_bayes(num, tab)
        elif self.tab_caption == 'Classifier':    
            self.form_classifier(num, tab)
        else:
            print 'Choice not found.'


    def form_root(self, num, tab):
        # tabFrame, row = self.set_inFrame(tab)
        tabFrame = tab
        row = 0
        row = self.entry_main_label(tabFrame, row, num)
        
        row = self.Entry(tabFrame, 'Root Fasta:', self.rootFasta_var, self.rootFasta, row,lenInputFrame10)
        row = self.Entry(tabFrame, 'Root Image:', self.rootImage_var, self.rootImage, row,lenInputFrame10)
        row = self.Entry(tabFrame, 'Root Files:', self.rootTable_var, self.rootTable, row,lenInputFrame10)
        row = self.Entry(tabFrame, 'Root Trees:', self.rootTree_var, self.rootTree, row,lenInputFrame10)
        row = self.Entry(tabFrame, 'Root Entropy Tables:', self.rootEntropy_var, self.rootEntropy, row,lenInputFrame10)

        row = self.break_line(tabFrame,  row)

        row = self.Entry(tabFrame, 'Mr.Bayes root:', self.mrBayes_path_var, self.mrBayes_path, row,lenInputFrame10)
        row = self.Entry(tabFrame, 'FigTree exec.:', self.figtree_path_var, self.figtree_path, row,lenInputFrame10)

        tk.Label(tabFrame, text="Standard",font=("Arial",10,"bold")).grid(row=row, column=1, sticky=tk.W)
        row += 1
        row = self.Entry(tabFrame, 'Filename min/alig:', self.mrBayes_fn01_var, self.mrBayes_fn01, row,lenInputFrame10)
        row = self.Entry(tabFrame, 'Filename min/cons:', self.mrBayes_fn02_var, self.mrBayes_fn02, row,lenInputFrame10)
        row = self.Entry(tabFrame, 'Filename max/alig:', self.mrBayes_fn03_var, self.mrBayes_fn03, row,lenInputFrame10)
        row = self.Entry(tabFrame, 'Filename max/cons:', self.mrBayes_fn04_var, self.mrBayes_fn04, row,lenInputFrame10)

        tk.Label(tabFrame, text="Covarion",font=("Arial",10,"bold")).grid(row=row, column=1, sticky=tk.W)
        row += 1
        row = self.Entry(tabFrame, 'Filename min/alig:', self.mrBayes_fn01_cov_var, self.mrBayes_fn01_cov, row,lenInputFrame10)
        row = self.Entry(tabFrame, 'Filename min/cons:', self.mrBayes_fn02_cov_var, self.mrBayes_fn02_cov, row,lenInputFrame10)
        row = self.Entry(tabFrame, 'Filename max/alig:', self.mrBayes_fn03_cov_var, self.mrBayes_fn03_cov, row,lenInputFrame10)
        row = self.Entry(tabFrame, 'Filename max/cons:', self.mrBayes_fn04_cov_var, self.mrBayes_fn04_cov, row,lenInputFrame10)


        return row

   
    def updtcblist(self):
        self.searchBox['values'] = self.search_list        

    def entry_main_label(self, inFrame, row, num):
        row = self.break_line( inFrame,  row)
        tk.Label(inFrame, text=self.text_list[num][1],font=("Arial",10,"bold")).grid(row=row, column=1, sticky=tk.W)
        row += 1
        row = self.break_line(inFrame, row)
        
        # print num, self.text_list[num][1]

        if num == 0:    
            f1 = tk.Frame(inFrame)
    
            self.searchBox = ttk.Combobox(f1, textvariable=self.search_org_gene_var, width=40, state='readonly', postcommand = self.updtcblist)
            self.searchBox['values'] = self.search_list
            
            self.searchBox.bind('<<ComboboxSelected>>', self.cmd_search_list)
            self.searchBox.set(self.first)
                
            self.searchBox.grid(column=0, row=0, sticky=tk.W)
        
            self.Button(f1, "New Search", self.clear_screen, row=0, col=1)
            self.Button(f1, "Delete", self.delete_org_gene, row=0, col=2)  
                
            f1.grid(row=row, column=1, sticky=tk.W)
            # self.row_cmd_search_list = row               
        elif num == 2:
            
            f1 = tk.Frame(inFrame)
            
            self.Button(f1, "SeaView: mincut", self.call_seaview_min, row=0, col=1)   
            self.Button(f1, "SeaView: maxmer", self.call_seaview_max, row=0, col=2)   
            self.Button(f1, "SeaView: species", self.call_seaview_species, row=0, col=3)
            
            f1.grid(row=row, column=1, sticky=tk.W)

        row += 1

                        
        return row
    
    def call_seaview_min(self):
        filename = self.rootFasta + self.aligned_filename_fasta_var.get().replace("(max/min)","mincut")
        self.call_external_program(self.seaview, filename)
        
    def call_seaview_max(self):
        filename = self.rootFasta + self.aligned_filename_fasta_var.get().replace("(max/min)","maxmer")
        self.call_external_program(self.seaview, filename)


    def call_seaview_species(self):
        if self.species == '':
            tkMessageBox.showinfo('Warning', 'Define a species, click in the species list')
            return
            
        filename = self.rootFasta + '%s_%s_%s_%s_%iL_aligned.fasta' %\
                 (self.organism, self.species, self.seqType, self.gene_title, self.cutoffLength)
                         
        self.call_external_program(self.seaview, filename) 
 
    def break_line(self, inFrame,  row):
        tk.Label(inFrame, text='').grid(row=row)
        row += 1
        return row
        
        
    def Choose_MinMax(self, inFrame, row):
        f1 = tk.Frame(inFrame, width=60, borderwidth=3, relief=tk.RAISED)
        self.seqType_var.set(self.seqType) # initialize
    
        tk.Radiobutton(f1, text="min", variable=self.minmax_var, value="mincut", indicatoron=0, command = self.on_change_fields).grid(row=0, column=2)
        tk.Radiobutton(f1, text="max", variable=self.minmax_var, value="maxmer", indicatoron=0, command = self.on_change_fields).grid(row=0, column=3)
        f1.grid(row=row, column=0, sticky=tk.W)

        row += 1
        return row
    
    def on_change_each_all(self):
        if self.each_all_var.get() == "all":
            if tkMessageBox.askyesno("Parameters", "May I change the parameters do default?"):
                self.showGraph = False
                self.saveGraph = True
                self.saveData = True
                
                self.showGraph_var.set(self.showGraph)
                self.saveGraph_var.set(self.saveGraph)
                self.saveData_var.set(self.saveData)
        
        
    def Choose_MinMax_CalcAll_option(self, inFrame, row):
        f1 = tk.Frame(inFrame, width=60, borderwidth=3, relief=tk.RAISED)
        self.seqType_var.set(self.seqType) # initialize
    
        tk.Radiobutton(f1, text="min", variable=self.minmax_var, value="mincut", indicatoron=0, command = self.on_change_fields).grid(row=0, column=2)
        tk.Radiobutton(f1, text="max", variable=self.minmax_var, value="maxmer", indicatoron=0, command = self.on_change_fields).grid(row=0, column=3)

        f2 = tk.Frame(inFrame, width=60, borderwidth=3, relief=tk.RAISED)
        tk.Radiobutton(f2, text="calc each", variable=self.each_all_var, value="each", indicatoron=0, command = self.on_change_each_all).grid(row=0, column=4)
        tk.Radiobutton(f2, text="calc all", variable=self.each_all_var, value="all", indicatoron=0, command = self.on_change_each_all).grid(row=0, column=5)

        f1.grid(row=row, column=0, sticky=tk.W)
        f2.grid(row=row, column=1, sticky=tk.W)

        row += 1
        return row
    
                    
    def Entry_Button(self, inFrame, label, var, field, row, width=50, readonly=True, cmdLabel="", cmd=None):
        tk.Label(inFrame, text=label).grid(row=row, column=0, sticky=tk.W)
        
        if not readonly:
            tk.Entry(inFrame, textvariable=var, width=width, state='normal').grid(row=row, column=1, sticky=tk.W)
        else:
            tk.Entry(inFrame, textvariable=var, width=width, state='disabled').grid(row=row, column=1, sticky=tk.W)

        var.set(field)
        
        if cmd:
            tk.Button(inFrame, text=cmdLabel, command=cmd).grid(row=row, column=2, sticky=tk.W)

        
        row += 1
        return row
    
    
    def Entry(self, inFrame, label, var, field, row, width=10, readonly=False):
        tk.Label(inFrame, text=label).grid(row=row, column=0, sticky=tk.W)
        
        if not readonly:
            tk.Entry(inFrame, textvariable=var, width=width, state='normal').grid(row=row, column=1, sticky=tk.W)
        else:
            tk.Entry(inFrame, textvariable=var, width=width, state='disabled').grid(row=row, column=1, sticky=tk.W)
        
        var.set(field)
        tk.Label(inFrame, text='', width=secondFieldWidth).grid(row=row, column=2)
        row += 1
        return row

    def Entry2lables(self, inFrame, label, label2, var, field, row, width):
        tk.Label(inFrame, text=label).grid(row=row, column=0, sticky=tk.W)
        
        f1 = tk.Frame(inFrame)
        tk.Entry(f1, textvariable=var, width=width).grid(row=row, column=0, sticky=tk.W)
        tk.Label(f1, text=label2).grid(row=row, column=1, sticky=tk.W)
        f1.grid(row=row, column=1, sticky=tk.W)

        
        tk.Label(inFrame, text='', width=secondFieldWidth).grid(row=row, column=2)
        
        var.set(field)
        row += 1
        return row

    def Entry2(self, inFrame, label, var, field, row, col, width):
        row -= 1
        tk.Label(inFrame, text=label).grid(row=row, column=col, sticky=tk.W)
        tk.Entry(inFrame, textvariable=var, width=width).grid(row=row, column=col+1, sticky=tk.W)
        var.set(field)
        tk.Label(inFrame, text='', width=secondFieldWidth).grid(row=row, column=2)
        row += 1
        return row
    
    def EntryRow(self, inFrame, label, var, field, row, col, width):
        tk.Label(inFrame, text=label).grid(row=row, column=col, sticky=tk.W)
        tk.Entry(inFrame, textvariable=var, width=width).grid(row=row, column=col+1, sticky=tk.W)
        var.set(field)
        tk.Label(inFrame, text='', width=secondFieldWidth).grid(row=row, column=2)
        return

    def EntryFrame(self, inFrame, label, var, field, row, col,  width=10, increment=False):
        tk.Label(inFrame, text=label).grid(row=row, column=col, sticky=tk.W)
        tk.Entry(inFrame, textvariable=var, width=width).grid(row=row, column=col+1, sticky=tk.W)
        var.set(field)
        if increment:
            row += 1
        return row
        
    def Check(self, inFrame, label, labelCheck, var, field, row):
        tk.Label(inFrame, text=label).grid(row=row, column=0, sticky=tk.W)
        tk.Checkbutton(inFrame, variable=var, text=labelCheck, command=None, onvalue=True, offvalue=False).grid(row=row, column=1, sticky=tk.W)
        var.set(field)
        row += 1
        return row         
        
    def ComboEntry(self, inFrame, label, lista, var, field, row, col=0, width=10):
        tk.Label(inFrame, text=label).grid(row=row, column=col, sticky=tk.W)

        comboList = ttk.Combobox(inFrame, textvariable=var, state='readonly', width=width)
        comboList['values'] = lista
        
        comboList.set(field)
        # comboList.current(0)
        comboList.grid(column=col+1, row=row, sticky=tk.W)
        row += 1
        return row
    
    
    def Checkbutton_Entry(self, frame, label, check_var, check_field, command,var, field, row, width):
        tk.Checkbutton(frame, variable=check_var, text=label, command=command, onvalue=1, offvalue=0).grid(row=row, column=0, sticky=tk.W)
        check_var.set(check_field)

        tk.Entry(frame, textvariable=var, width=width).grid(row=row, column=1, sticky=tk.W)
        var.set(field)
        tk.Label(frame, text='', width=secondFieldWidth).grid(row=row, column=2)
        row += 1
        return row

    def Checkbutton(self, frame, label, var, field, col, command=None):
        if command:
            tk.Checkbutton(frame, text=label, variable=var,command=command).grid(row=0, column=col)
        else:
            tk.Checkbutton(frame, text=label, variable=var).grid(row=0, column=col)
        var.set(field) 
        

    def Button(self, frame, label, command, row, col):
        tk.Button(frame, text=label, command=command).grid(row=row, column=col, sticky=tk.W)
        row += 1
        return row


    def Button_Vert_Horiz(self, frame, label, row):
        tk.Label(frame, text=label).grid(row=row, column=0, sticky=tk.W)

        f1 = tk.Frame(frame, width=72, borderwidth=3, relief=tk.RAISED)
        self.vert_horiz_var.set(self.vert_horiz) # initialize
            
        tk.Radiobutton(f1, text="Vertical Entropy", variable=self.vert_horiz_var, value="VSH", indicatoron=0, command = self.rb_vert_horiz_selected).grid(row=0, column=0)
        tk.Radiobutton(f1, text="Vertical MI", variable=self.vert_horiz_var, value="VMI", indicatoron=0, command = self.rb_vert_horiz_selected).grid(row=0, column=1)
        tk.Radiobutton(f1, text="Horizontal MI", variable=self.vert_horiz_var, value="HMI", indicatoron=0, command = self.rb_vert_horiz_selected).grid(row=0, column=2)

        tk.Label(f1, text="      ").grid(row=0, column=3, sticky=tk.W)

        tk.Radiobutton(f1, text="standard", variable=self.rand_method_var, value="standard", indicatoron=0, command = self.on_change_fields).grid(row=0, column=4)
        tk.Radiobutton(f1, text="shuffling", variable=self.rand_method_var, value="shuffle", indicatoron=0, command = self.on_change_fields).grid(row=0, column=5)
        tk.Radiobutton(f1, text="randomize", variable=self.rand_method_var, value="random", indicatoron=0, command = self.on_change_fields).grid(row=0, column=6)

        f1.grid(row=row, column=1, sticky=tk.W)
        row += 1
        return row

    
    def Button_Cluster_Methods(self, frame, label, row):
        tk.Label(frame, text=label).grid(row=row, column=0, sticky=tk.W)

        f1 = tk.Frame(frame, width=60, borderwidth=3, relief=tk.RAISED)
        self.cluster_method_var.set(self.cluster_method) # initialize
            
        tk.Radiobutton(f1, text="Complete", variable=self.cluster_method_var, value="complete", indicatoron=0, command = self.rb_vert_horiz_selected).grid(row=0, column=0)
        tk.Radiobutton(f1, text="Single", variable=self.cluster_method_var, value="single", indicatoron=0, command = self.rb_vert_horiz_selected).grid(row=0, column=1)
        tk.Radiobutton(f1, text="WPGMA", variable=self.cluster_method_var, value="weighted", indicatoron=0, command = self.rb_vert_horiz_selected).grid(row=0, column=2)
        tk.Radiobutton(f1, text="Centroid", variable=self.cluster_method_var, value="centroid", indicatoron=0, command = self.rb_vert_horiz_selected).grid(row=0, column=3)

        f1.grid(row=row, column=1, sticky=tk.W)
        row += 1
        return row

    def three_options_yes_no(self, frame, lable, row, column, var, default, command=None):
        if column == 5:
            tk.Label(frame, text='').grid(row=row, column=4)
       
        tk.Label(frame, text=lable).grid(row=row, column=column, sticky=tk.W)

        var.set(default) # initialize
        tk.Radiobutton(frame, text="no", variable=var, value="no", indicatoron=0, command = command).grid(row=row, column=column+1)
        tk.Radiobutton(frame, text="yes", variable=var, value="yes", indicatoron=0, command = command).grid(row=row, column=column+2)
        tk.Radiobutton(frame, text="---", variable=var, value="---", indicatoron=0, command = command).grid(row=row, column=column+3)

        
    def form_NCBI(self, num, tab):
        # inFrame, row = self.set_inFrame(tab)
        inFrame = tab
        row = 0
        
        row = self.entry_main_label(inFrame, row, num)
        
        row = self.Entry(inFrame, 'Organism:', self.organism_var, self.organism, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Gene list:', self.geneList_var, self.geneList, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Gene:', self.gene_var, self.gene, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Species list:', self.specieList_var, self.specieList, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Title:', self.title_var, self.title, row,lenInputFrame10)
        
        lista = ['Nucleotide','Protein']
        ''' row = self.Entry(inFrame, 'Database:', self.db_var, self., row,lenInputFrame10)
            ComboEntry(self, inFrame, label, lista, var, field, row) '''
        row = self.ComboEntry(inFrame, 'Database:', lista, self.db_var, self.db, row, col=0, width=10)

        row = self.break_line( inFrame,  row)
        
        '''----------------------------------------------------'''

        tk.Label(inFrame, text='Title filters:').grid(row=row, column=0, sticky=tk.W)
        f1 = tk.Frame(inFrame, width=60, borderwidth=3, relief=tk.RAISED)
        f1.grid(row=row, column=1, sticky=tk.W)
        
        self.three_options_yes_no(f1, 'CDS:', row=row, column=0, var=self.cds_var, default=self.cds, command=self.rb_seq_yes_no_selected)
        self.three_options_yes_no(f1, 'mRNA  ', row=row, column=5, var=self.mrna_var, default=self.mrna)
        self.three_options_yes_no(f1, 'Partial CDS:', row=row, column=10, var=self.partialCDS_var, default=self.partialCDS)

        '''----------------------------------------------------'''
        row += 1
        self.three_options_yes_no(f1, 'Complete Sequence:', row=row, column=0, var=self.completeSeq_var, default=self.completeSeq, command=self.rb_seq_yes_no_selected)
        self.three_options_yes_no(f1, 'Complete Genome:  ', row=row, column=5, var=self.completeGen_var, default=self.shotGun)
        self.three_options_yes_no(f1, 'Shotgun:        ', row=row, column=10, var=self.shotGun_var, default=self.shotGun)
        
        '''----------------------------------------------------'''
        row += 1
        self.three_options_yes_no(f1, 'Chromosome:     ', row=row, column=0, var=self.chromosome_var, default=self.chromosome)
        self.three_options_yes_no(f1, 'Contig:           ', row=row, column=5, var=self.contig_var, default=self.contig)
        self.three_options_yes_no(f1, 'Super-Contig:     ', row=row, column=10, var=self.superContig_var, default=self.superContig)

        '''----------------------------------------------------'''
        row = self.break_line(inFrame, row)
        
        row = self.Entry(inFrame, 'Log Gbk Filename:', self.logGbkFilename_var, self.logGbkFilename, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Output:', self.output_filename_gbk_var, self.output_filename_gbk, row,lenInputFrame10)
        row = self.Entry(inFrame, 'retmax:', self.retmax_var, self.retmax, row,10)
        row = self.Entry(inFrame, 'e-mail:', self.email_var, self.email, row,lenInputFrame10)
        #row = self.break_line( inFrame,  row)
        row = self.Entry(inFrame, 'Message:', self.msg_obs_var, self.msg_obs, row,lenInputFrame10)

        # row = self.break_line( inFrame,  row)

            
    def clear_screen(self):
        self.new_search = True
        
        self.organism = ""
        self.gene = ""
        self.title = ""
        self.geneList = ""
        self.specieList = ""
        
        self.gene_title = ""
        
        self.organism2 = ""
        self.gene2 = ""
        self.title2 = ""

        self.dicParams = {}
        
        self.organism_var.set(self.organism)
        self.gene_var.set(self.gene)
        self.title_var.set(self.title)
        self.geneList_var.set(self.geneList)
        self.specieList_var.set(self.specieList)
        
        self.speciesListbox.delete(0, tk.END) 
        
        self.clear_text_area()
        
        self.showmsg_obs("Set organism, gene, if necessary gene list and/or title, and run get GBK from NCBI")
        
                
    def form_GBK_Fasta(self, num, tab):
        # inFrame, row = self.set_inFrame(tab)
        inFrame = tab
        row = 0
        
        row = self.entry_main_label( inFrame,  row, num)
        
        row = self.Entry(inFrame, 'Organism:', self.organism_var, self.organism, row,lenInputFrame10,readonly=True)
        #row = self.Entry(inFrame, 'Gene list:', self.geneList_var, self.geneList, row,100)
        row = self.Entry(inFrame, 'Gene:', self.gene_var, self.gene, row,lenInputFrame10,readonly=True)
        row = self.Entry(inFrame, 'Title:', self.title_var, self.title, row,lenInputFrame10,readonly=True)

        row = self.break_line(inFrame,  row)
        row = self.Entry(inFrame, 'Cutoff Length:', self.cutoffLength_var, self.cutoffLength, row,10)
        row = self.Entry(inFrame, 'Maxi length seq:', self.maxiSeq_var, self.maxiSeq, row,10)
        row = self.Entry(inFrame, 'Mini length seq:', self.minSeq_var, self.minSeq, row,10)

        row = self.break_line(inFrame,  row)
        f2 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        f2.grid(row=row, column=1, sticky=tk.W)
        row += 1
        
        self.Checkbutton( f2, "Stop", self.stop_var, self.stop, col=0)        
        self.Checkbutton( f2, "Show Message", self.showmessage_var, self.showmessage, col=1)        
        
        row = self.break_line( inFrame,  row)

        row = self.Entry(inFrame, 'Input', self.output_filename_gbk_var, self.output_filename_gbk, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Output', self.output_filename_fasta_var, self.output_filename_fasta, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Message:', self.msg_obs_var, self.msg_obs, row,lenInputFrame10)

        row = self.break_line( inFrame,  row)


    def form_aligment(self, num, tab):
        # inFrame, row = self.set_inFrame(tab)
        inFrame = tab
        row = 0

        row = self.entry_main_label( inFrame,  row, num)
        
        row = self.Entry(inFrame, 'Organism:', self.organism_var, self.organism, row,lenInputFrame10,readonly=True)
        row = self.Entry(inFrame, 'Gene:', self.gene_var, self.gene, row,lenInputFrame10,readonly=True)
        row = self.Entry(inFrame, 'Title:', self.title_var, self.title, row,lenInputFrame10,readonly=True)
        row = self.Entry(inFrame, 'Species:', self.species_var, self.species, row,lenInputFrame10)

        row = self.break_line(inFrame,  row)
        row = self.Entry(inFrame ,'Cutoff #Seqs',   self.cutoffNumSeq_var, self.cutoffNumSeq, row,5)
        row = self.Entry(inFrame, 'Cutoff Length:', self.cutoffLength_var, self.cutoffLength, row,5)

        '''
          Muscle (MUltiple Sequence Comparison by Log-Expectation) - http://www.drive5.com/muscle/
          Clustal  - www.clustal.org/omega/
        '''
        #row = self.Entry(inFrame, 'Vertical Cutoff (%):', self.horizCutoff_var, self.horizCutoff, row,10)
        listaAlignment = ['Muscle']  # ,'Clustal'
        row = self.ComboEntry(inFrame, 'Method:', listaAlignment, self.alignment_var, self.alignment, row, col=0, width=10)

        row += 1
        # row = self.break_line(inFrame, row)

        row = self.Entry2lables(inFrame, 'Max Horizontal gaps:', "%", self.horizCutoff_var, self.horizCutoff, row,5)
        row = self.Entry2lables(inFrame, 'Max Vertical gaps1:', "% low value creates mincut", self.minVertCutoff_var, self.minVertCutoff, row,5)
        row = self.Entry2lables(inFrame, 'Max Vertical gaps2:', "% high value creates maxmer", self.maxVertCutoff_var, self.maxVertCutoff, row,5)

        '''-----------------------------------------'''
        f2 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        f2.grid(row=row, column=1, sticky=tk.W)
        
        self.Checkbutton(f2, "Realign each species", self.realing_each_species_var, self.realing_each_species, col=1)        
        self.Checkbutton(f2, "Realign all together", self.realing_all_var, self.realing_all, col=2)        
        row += 1
        
        row = self.break_line(inFrame, row)

        row = self.Entry(inFrame, 'Input', self.output_filename_fasta_var, self.output_filename_fasta, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Output', self.aligned_filename_fasta_var, self.aligned_filename_fasta, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Message:', self.msg_obs_var, self.msg_obs, row,lenInputFrame10)

        
    def form_Purging(self, num, tab):
        # inFrame, row = self.set_inFrame(tab)
        inFrame = tab
        row = 0

        row = self.entry_main_label(inFrame, row, num)
                    
        row = self.Entry(inFrame, 'Organism:', self.organism_var, self.organism, row,lenInputFrame70,readonly=True)
        row = self.Entry(inFrame, 'Gene:', self.gene_var, self.gene, row,lenInputFrame70,readonly=True)
        row = self.Entry(inFrame, 'Title:', self.title_var, self.title, row,lenInputFrame10,readonly=True)
        row = self.Entry(inFrame ,'Cutoff #Seqs',   self.cutoffNumSeq_var, self.cutoffNumSeq, row,5,readonly=True)
        row = self.Entry(inFrame, 'Cutoff Length:', self.cutoffLength_var, self.cutoffLength, row,5,readonly=True)
  
        row = self.break_line(inFrame, row)
 
        row = self.entry_seqType(inFrame, row, False, False) 
             
        row = self.break_line(inFrame, row)

        f2 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        f2.grid(row=row, column=1, sticky=tk.W)
        
        self.Checkbutton(f2, "Confirm gene", self.confirm_gene_var, self.confirm_gene, col=1)        
        self.Checkbutton(f2, "Show Message", self.showmessage_var, self.showmessage, col=2)        
        row += 1
             
        row = self.break_line(inFrame, row)
        
        row = self.Entry(inFrame, 'List bad terms', self.listBadTerms_var, self.listBadTerms, row,lenInputFrame10)
             
        row = self.break_line(inFrame, row)

        row = self.Entry(inFrame, 'Input', self.aligned_filename_fasta_var, self.aligned_filename_fasta, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Output', self.output_filename_purged_var, self.output_filename_purged, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Message:', self.msg_obs_var, self.msg_obs, row, lenInputFrame10)

        row = self.break_line(inFrame,  row)

        ''' rever !!!
        row = self.Button(inFrame, "Merge all aligned, creates maxmer (mincut is manual)", command=self.run_merge_all, row=row, col=1)
        '''
        # row = self.Button(inFrame, "Hierarchical Cluster", command=self.run_hierarchical_cluster, row=row, col=1)        

    def entry_seqType(self, inFrame, row, begin=True, end=True):
        if begin:        
            row = self.break_line(inFrame, row)

        f1 = tk.Frame(inFrame, width=60, borderwidth=3, relief=tk.RAISED)
        self.seqType_var.set(self.seqType) # initialize
    
        tk.Radiobutton(f1, text="Gene", variable=self.seqType_var, value="Gene", indicatoron=0, command = self.rb_seqType_selected).grid(row=0, column=0)
        tk.Radiobutton(f1, text="CDS", variable=self.seqType_var, value="CDS", indicatoron=0, command = self.rb_seqType_selected).grid(row=0, column=1)
        tk.Radiobutton(f1, text="Exon", variable=self.seqType_var, value="Exon", indicatoron=0, command = self.rb_seqType_selected).grid(row=0, column=2)

        f1.grid(row=row, column=1, sticky=tk.W)
        row += 1

        if end:        
            row = self.break_line(inFrame, row)

        return row
    
    def entry_seqType_random(self, inFrame, row, begin=True, end=True):
        if begin:        
            row = self.break_line(inFrame, row)

        f1 = tk.Frame(inFrame, width=60, borderwidth=3, relief=tk.RAISED)
        self.seqType_var.set(self.seqType) # initialize
    
        tk.Radiobutton(f1, text="Gene", variable=self.seqType_var, value="Gene", indicatoron=0, command = self.rb_seqType_selected).grid(row=0, column=0)
        tk.Radiobutton(f1, text="CDS", variable=self.seqType_var, value="CDS", indicatoron=0, command = self.rb_seqType_selected).grid(row=0, column=1)
        tk.Radiobutton(f1, text="Exon", variable=self.seqType_var, value="Exon", indicatoron=0, command = self.rb_seqType_selected).grid(row=0, column=2)


        tk.Label(f1, text="      ").grid(row=0, column=3, sticky=tk.W)

        tk.Radiobutton(f1, text="standard", variable=self.rand_method_var, value="standard", indicatoron=0, command = self.on_change_fields).grid(row=0, column=4)
        tk.Radiobutton(f1, text="shuffling", variable=self.rand_method_var, value="shuffle", indicatoron=0, command = self.on_change_fields).grid(row=0, column=5)
        tk.Radiobutton(f1, text="randomize", variable=self.rand_method_var, value="random", indicatoron=0, command = self.on_change_fields).grid(row=0, column=6)

        f1.grid(row=row, column=1, sticky=tk.W)
        row += 1

        if end:        
            row = self.break_line(inFrame, row)

        return row

    def onKeyPress(self, event):
        print 'Got key press:', event.char, str(event)
        
    def callback(self, event):
        print "clicked at", event.x, event.y
       
    def on_change_fields(self, event=None):
        self.get_params()

        if self.rand_method_var.get() == "shuffle":
            self.stri_random = "_shuffle"
            self.label_random = "shuffle"
            
        elif self.rand_method_var.get() == "random":
            self.stri_random = "_random"
            self.label_random = "random"
            
        else:
            self.stri_random = ""
            self.label_random = ""
        
        self.species_var.set(self.species)

        self.set_names()                
        
        self.output_filename_gbk_var.set(self.output_filename_gbk)
        self.output_filename_fasta_var.set(self.output_filename_fasta)
        self.aligned_filename_fasta_var.set(self.aligned_filename_fasta)
        self.input_filename_align_var.set(self.input_filename_align)
        self.output_filename_purged_var.set(self.output_filename_purged)
        self.output_filename_consensus_var.set(self.output_filename_consensus)
        self.output_hmi_var.set(self.output_hmi)
        self.output_vmi_var.set(self.output_vmi)
        self.logGbkFilename_var.set(self.logGbkFilename)
        
        self.input_jsd_var.set(self.input_jsd)
        self.cluster_input_filename_var.set(self.cluster_input_filename)


    def form_Consensus(self, num, tab):
        # inFrame, row = self.set_inFrame(tab)
        inFrame = tab
        row = 0
        
        row = self.entry_main_label(inFrame, row, num)

        row = self.Entry(inFrame, 'Organism:', self.organism_var, self.organism, row,lenInputFrame70,readonly=True)
        row = self.Entry(inFrame, 'Gene:', self.gene_var, self.gene, row,lenInputFrame70,readonly=True)
        row = self.Entry(inFrame, 'Title:', self.title_var, self.title, row,lenInputFrame10,readonly=True)

        row = self.break_line(inFrame, row)

        row = self.Entry(inFrame ,'Cutoff #Seqs',   self.cutoffNumSeq_var, self.cutoffNumSeq, row,5,readonly=True)
        row = self.Entry(inFrame, 'Cutoff Length:', self.cutoffLength_var, self.cutoffLength, row,5,readonly=True)

        row = self.entry_seqType(inFrame, row, False, False) 

        f2 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        f2.grid(row=row, column=1, sticky=tk.W)
        tk.Checkbutton(f2, text="Start and Stop Codon", variable=self.start_stop_codon_var).grid(row=0, column=0)
        self.start_stop_codon_var.set(self.start_stop_codon)   
        
        tk.Checkbutton(f2, text="Show Message", variable=self.showmessage_var).grid(row=0, column=1)
        row += 1

        row = self.break_line(inFrame, row)  

        row = self.Entry(inFrame, 'Input', self.output_filename_purged_var, self.output_filename_purged, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Output', self.output_filename_consensus_var, self.output_filename_consensus, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Message:', self.msg_obs_var, self.msg_obs, row, lenInputFrame10)

        row = self.break_line(inFrame, row)


    def form_HMI(self, num, tab):
        # inFrame, row = self.set_inFrame(tab)
        inFrame = tab
        row = 0
        
        row = self.entry_main_label( inFrame,  row, num)
        
        row = self.Entry(inFrame, 'Organism:', self.organism_var, self.organism, row,lenInputFrame70,readonly=True)
        row = self.Entry(inFrame, 'Gene:', self.gene_var, self.gene, row,lenInputFrame70,readonly=True)
        row = self.Entry(inFrame, 'Title:', self.title_var, self.title, row,lenInputFrame10,readonly=True)
        row = self.Choose_MinMax_CalcAll_option(inFrame, row)

        row = self.break_line( inFrame,  row)
        
        row = self.Entry(inFrame, 'Num of letters:', self.numOfLetters_var, self.numOfLetters, row,10)
        row = self.Entry(inFrame, 'Frame:', self.frame_var, self.frame, row,10)
        row = self.Entry(inFrame, 'Offset:', self.offset_var, self.offset, row,10)
        row = self.entry_seqType_random(inFrame, row, False, False) 

        f1a = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        self.Checkbutton(f1a, "Save data", self.saveData_var, self.saveData, col=0)        
        self.Checkbutton(f1a, "See image", self.showGraph_var, self.showGraph, col=1)        
        self.Checkbutton(f1a, "Save image", self.saveGraph_var, self.saveGraph, col=2)        
        #self.Checkbutton(f1a, "Log plot", self.isLog_var, self.isLog, col=3)        
        self.Checkbutton(f1a, "Recalculate", self.de_novo_var, self.de_novo, col=4)        
        self.Checkbutton(f1a, "Normalized", self.norm_var, self.norm, col=5)    
        f1a.grid(row=row, column=1, sticky=tk.W)
        row += 1

        f1b = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        self.Checkbutton( f1b, "bias correction", self.withCorrection_var, self.withCorrection, col=0, command=self.on_change_fields)        
        self.Checkbutton( f1b, "mnat", self.mnat_var, self.mnat, col=1)
        f1b.grid(row=row, column=1, sticky=tk.W)
        row += 1  

        f2 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        f2.grid(row=row, column=1, sticky=tk.W)
        row = self.ComboEntry(f2, 'Image: ', self.listImage, self.imgType_var, self.imgType, row)
        row = self.Entry2(f2, 'DPI: ', self.dpi_var, self.dpi, row, col=3, width=5)
        row += 1  
        
        row = self.break_line( inFrame,  row)

        row = self.Entry(inFrame, 'Input', self.output_filename_consensus_var, self.output_filename_consensus, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Output', self.output_hmi_var, self.output_hmi, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Message:', self.msg_obs_var, self.msg_obs, row, lenInputFrame10)

        row = self.break_line( inFrame,  row)
        
        
    def form_VMI(self, num, tab):
        # inFrame, row = self.set_inFrame(tab)
        inFrame = tab
        row = 0
        
        row = self.entry_main_label( inFrame,  row, num)
        
        row = self.Entry(inFrame, 'Organism:', self.organism_var, self.organism, row,lenInputFrame70,readonly=True)
        row = self.Entry(inFrame, 'Gene:', self.gene_var, self.gene, row,lenInputFrame70,readonly=True)
        row = self.Entry(inFrame, 'Title:', self.title_var, self.title, row,lenInputFrame10,readonly=True)
        row = self.Choose_MinMax_CalcAll_option(inFrame, row)
        
        row = self.break_line(inFrame, row)

        row = self.Entry(inFrame, 'Num of letters:', self.numOfLetters_var, self.numOfLetters, row,10)
        row = self.entry_seqType_random(inFrame, row, False, False) 
        
        f2 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        f2.grid(row=row, column=1, sticky=tk.W)
        self.Checkbutton( f2, "Save data", self.saveData_var, self.saveData,  col=0)        
        self.Checkbutton( f2, "See image", self.showGraph_var, self.showGraph, col=1)        
        self.Checkbutton( f2, "Save image", self.saveGraph_var, self.saveGraph, col=2)        
        #self.Checkbutton( f2, "Log plot", self.isLog_var, self.isLog, col=3)        
        self.Checkbutton( f2, "Recalculate", self.de_novo_var, self.de_novo, col=4)    
        self.Checkbutton( f2, "Normalized", self.norm_var, self.norm, col=5)    
        row += 1

        f1 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        f1.grid(row=row, column=1, sticky=tk.W)
        self.Checkbutton( f1, "bias correction", self.withCorrection_var, self.withCorrection, col=0, command=self.on_change_fields)        
        self.Checkbutton( f1, "mnat", self.mnat_var, self.mnat, col=1)
        row += 1  
        
        f2 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        f2.grid(row=row, column=1, sticky=tk.W)
        row = self.ComboEntry(f2, 'Image: ', self.listImage, self.imgType_var, self.imgType, row)
        row = self.Entry2(f2, 'DPI: ', self.dpi_var, self.dpi, row, col=3, width=5)

        row = self.break_line( inFrame,  row)

        f3 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        tk.Label(f3, text="Heatmap 2D").grid(row=0, column=0, sticky=tk.W)
        tk.Checkbutton(f3, text="apply", variable=self.heatmap_ceil_var).grid(row=0, column=1)        
        self.EntryFrame(f3, 'ceil:', self.heatmap_ceil_value_var, self.heatmap_ceil_value, row=0, col=2, width=6)
        self.Checkbutton(f3, "3D", self.is3D_var, self.is3D, col=4)        
        self.ComboEntry(f3, 'color scheme: ', self.listHeatMapColors, self.heatmap_color_var, self.heatmap_color, row=0, col=5, width=8)
        f3.grid(row=row, column=1, sticky=tk.W)
        row += 1


        row = self.break_line( inFrame,  row)

        row = self.Entry(inFrame, 'Input', self.output_filename_consensus_var, self.output_filename_consensus, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Output', self.output_vmi_var, self.output_vmi, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Message:', self.msg_obs_var, self.msg_obs, row, lenInputFrame10)

        row = self.break_line( inFrame,  row)
        
           
    def form_JSD(self, num, tab):
        # inFrame, row = self.set_inFrame(tab)
        inFrame = tab
        row = 0
        
        row = self.entry_main_label( inFrame,  row, num)
        
        row = self.Entry(inFrame, 'Organism:', self.organism_var, self.organism, row,lenInputFrame70,readonly=True)
        row = self.Entry(inFrame, 'Gene:', self.gene_var, self.gene, row,lenInputFrame70,readonly=True)
        row = self.Entry(inFrame, 'Title:', self.title_var, self.title, row,lenInputFrame10,readonly=True)
        row = self.Choose_MinMax_CalcAll_option(inFrame, row)
        
        row = self.Entry(inFrame, 'Num of letters:', self.numOfLetters_var, self.numOfLetters, row,10)
        
        row = self.entry_seqType(inFrame, row, False, False) 
        row = self.Entry(inFrame, '#SE:', self.numOfSEs_var, self.numOfSEs, row,10)

        #row = self.break_line( inFrame,  row)

        row = self.Button_Vert_Horiz(inFrame, 'Analysis Method:', row)
        row = self.Entry(inFrame, 'Frame:', self.frame_var, self.frame, row,10)

        #row = self.break_line( inFrame,  row)
               
        f2 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        f2.grid(row=row, column=1, sticky=tk.W)
        self.Checkbutton( f2, "Save data", self.saveData_var, self.saveData, col=0)        
        self.Checkbutton( f2, "See image", self.showGraph_var, self.showGraph, col=1)        
        self.Checkbutton( f2, "Save image", self.saveGraph_var, self.saveGraph, col=2)        
        #self.Checkbutton( f2, "Log plot", self.isLog_var, self.isLog, col=3)        
        self.Checkbutton( f2, "Recalculate", self.de_novo_var, self.de_novo, col=4)        
        row += 1

        
        f3 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        f3.grid(row=row, column=1, sticky=tk.W)
        self.Checkbutton( f3, "bias correction", self.withCorrection_var, self.withCorrection, col=0, command=self.on_change_fields)        
        self.Checkbutton( f3, "mnat", self.mnat_var, self.mnat, col=1)
        row += 1  
               
        f4 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        f4.grid(row=row, column=1, sticky=tk.W)
        row = self.ComboEntry(f4, 'Image: ', self.listImage, self.imgType_var, self.imgType, row)
        row = self.Entry2(f4, 'DPI: ', self.dpi_var, self.dpi, row, col=3, width=5)
        
        row = self.break_line( inFrame,  row)

        row = self.Entry(inFrame, 'Input', self.input_jsd_var, self.input_jsd, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Output', self.cluster_input_filename_var, self.cluster_input_filename, row,lenInputFrame10)
                         
        row = self.Entry(inFrame, 'Message:', self.msg_obs_var, self.msg_obs, row, lenInputFrame10)

        row = self.break_line( inFrame,  row)
        
        
    def form_cluster(self, num, tab):
        # inFrame, row = self.set_inFrame(tab)
        inFrame = tab
        row = 0
        
        row = self.entry_main_label( inFrame,  row, num)

        row = self.Entry(inFrame, 'Organism:', self.organism_var, self.organism, row,lenInputFrame70,readonly=True)
        row = self.Entry(inFrame, 'Gene:', self.gene_var, self.gene, row,lenInputFrame70,readonly=True)
        row = self.Entry(inFrame, 'Title:', self.title_var, self.title, row,lenInputFrame10,readonly=True)
        row = self.Choose_MinMax_CalcAll_option(inFrame, row)
        
        row = self.Entry(inFrame, 'Num of letters:', self.numOfLetters_var, self.numOfLetters, row,10)
        row = self.entry_seqType(inFrame, row, False, False) 

        # row = self.break_line( inFrame,  row)

        row = self.Button_Vert_Horiz(inFrame, 'Analysis Method:', row)
        row = self.Entry(inFrame, 'Frame:', self.frame_var, self.frame, row,10)
        row = self.Button_Cluster_Methods(inFrame, 'Method:', row)
        row = self.Entry(inFrame, 'Leaf threshold:', self.colorThreshold_var, self.colorThreshold, row,10)



        f2 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        f2.grid(row=row, column=1, sticky=tk.W)
        self.Checkbutton( f2, "See image", self.showGraph_var, self.showGraph, col=0)        
        self.Checkbutton( f2, "Save image", self.saveGraph_var, self.saveGraph, col=1)        
        #self.Checkbutton( f2, "Log plot", self.isLog_var, self.isLog, col=2)        
        row += 1

        
        f1 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        f1.grid(row=row, column=1, sticky=tk.W)
        self.Checkbutton( f1, "bias correction", self.withCorrection_var, self.withCorrection, col=0, command=self.on_change_fields)        
        self.Checkbutton( f1, "mnat", self.mnat_var, self.mnat, col=1)
        row += 1  
               
        f4 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        f4.grid(row=row, column=1, sticky=tk.W)
        row = self.ComboEntry(f4, 'Image: ', self.listImage, self.imgType_var, self.imgType, row)
        row = self.Entry2(f4, 'DPI: ', self.dpi_var, self.dpi, row, col=3, width=5)
        
        row = self.break_line( inFrame,  row)

        row = self.Entry(inFrame, 'Input', self.cluster_input_filename_var, self.cluster_input_filename, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Message:', self.msg_obs_var, self.msg_obs, row, lenInputFrame10)

        #row = self.break_line( inFrame,  row)

    def form_entropy(self, num, tab):
        # inFrame, row = self.set_inFrame(tab)
        inFrame = tab
        row = 0
        row = self.entry_main_label( inFrame,  row, num)

        lista = ['Nucleotide','Protein']
        row = self.ComboEntry(inFrame, 'Database:', lista, self.db_var, self.db, row, col=0, width=10)

        row = self.Entry(inFrame, 'Num of letters:', self.numOfLetters_var, self.numOfLetters, row,10)
        row = self.Entry(inFrame, 'Num of experiments:', self.numOfExperiments_var, self.numOfExperiments, row,10)
        row = self.Entry(inFrame, 'Start at:', self.sartAt_var, self.sartAt, row,10)
        row = self.Entry(inFrame, 'Length simulation:', self.lengthSim_var, self.lengthSim, row,10)

        f1 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        row = self.ComboEntry(f1, 'Image: ', self.listImage, self.imgType_var, self.imgType, row)
        row = self.Entry2(f1, 'DPI: ', self.dpi_var, self.dpi, row, col=3, width=5)
        f1.grid(row=row, column=1, sticky=tk.W)

        row = self.break_line( inFrame,  row)

        f2 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        self.Checkbutton( f2, "Save data", self.saveData_var, self.saveData, col=0)        
        self.Checkbutton( f2, "See image", self.showGraph_var, self.showGraph, col=1)        
        self.Checkbutton( f2, "Save image", self.saveGraph_var, self.saveGraph, col=2)        
        f2.grid(row=row, column=1, sticky=tk.W)
        row += 1

        row = self.break_line( inFrame,  row)

        row = self.Entry(inFrame, 'Output', self.shannon_filename_var, self.shannon_filename, row,lenInputFrame10)
        row = self.Entry(inFrame, 'Message:', self.msg_obs_var, self.msg_obs, row, lenInputFrame10)

        row = self.break_line( inFrame,  row)

                    
    def form_mr_bayes(self, num, tab):
        # inFrame, row = self.set_inFrame(tab)
        inFrame = tab
        row = 0
        row = self.entry_main_label( inFrame,  row, num)
        # row = self.Entry_Button(inFrame, 'Species:', self.species_var, self.species, row, lenInputFrame10, readonly=True, cmdLabel="Tree", cmd=self.build_see_tree)
        # row = self.Entry(inFrame, 'Species:', self.species_var, self.species, row, 100)

        f1 = tk.Frame(inFrame, width=60, borderwidth=3, relief=tk.RAISED)
        self.seqType_var.set(self.seqType) # initialize
    
        tk.Radiobutton(f1, text="min", variable=self.minmax_var, value="mincut", indicatoron=0, command = self.on_change_fields).grid(row=0, column=1)
        tk.Radiobutton(f1, text="max", variable=self.minmax_var, value="maxmer", indicatoron=0, command = self.on_change_fields).grid(row=0, column=2)
        
        tk.Label(f1, text="   ").grid(row=0, column=3, sticky=tk.W)
        tk.Radiobutton(f1, text="aligned", variable=self.aligned_consensus_var, value="aligned", indicatoron=0, command = self.on_change_fields).grid(row=0, column=4)
        tk.Radiobutton(f1, text="consensus", variable=self.aligned_consensus_var, value="consensus", indicatoron=0, command = self.on_change_fields).grid(row=0, column=5)

        tk.Label(f1, text="   ").grid(row=0, column=6, sticky=tk.W)
        tk.Radiobutton(f1, text="standard", variable=self.standard_covarion_var, value="standard", indicatoron=0, command = self.on_change_fields).grid(row=0, column=7)
        tk.Radiobutton(f1, text="covarion", variable=self.standard_covarion_var, value="covarion", indicatoron=0, command = self.on_change_fields).grid(row=0, column=8)
        
        tk.Label(f1, text="   ").grid(row=0, column=9, sticky=tk.W)
        self.Checkbutton(f1, "Collapsed tree", self.colaps_tree_var, self.colaps_tree, col=10)        

        tk.Label(f1, text="   ").grid(row=0, column=11, sticky=tk.W)
        tk.Button(f1, text="Tree", command=self.build_see_tree).grid(row=0, column=12)
        tk.Button(f1, text="FigTree", command=self.build_fig_tree).grid(row=0, column=13)

        f1.grid(row=row, column=1, sticky=tk.W)
        row += 1


        fmb1 = tk.Frame(inFrame, width=60, borderwidth=3, relief=tk.RAISED)
        self.Checkbutton(fmb1, "pi(A)", self.piA_var, self.piA, col=0) 
        self.Checkbutton(fmb1, "pi(C)", self.piC_var, self.piC, col=1) 
        self.Checkbutton(fmb1, "pi(G)", self.piG_var, self.piG, col=2) 
        self.Checkbutton(fmb1, "pi(T)", self.piT_var, self.piT, col=3) 
        
        self.Checkbutton(fmb1, "rAC", self.rAC_var, self.rAC, col=5) 
        self.Checkbutton(fmb1, "rAG", self.rAG_var, self.rAG, col=6) 
        self.Checkbutton(fmb1, "rAT", self.rAT_var, self.rAT, col=7) 
        self.Checkbutton(fmb1, "rCG", self.rCG_var, self.rCG, col=8) 
        self.Checkbutton(fmb1, "rCT", self.rCT_var, self.rCT, col=9) 
        self.Checkbutton(fmb1, "rGT", self.rGT_var, self.rGT, col=10) 
        fmb1.grid(row=row, column=1, sticky=tk.W)
        row += 1

     
        fmb3a = tk.Frame(inFrame, width=60, borderwidth=3, relief=tk.RAISED)
        self.Checkbutton(fmb3a, "LnL", self.LnL_var, self.LnL, col=0) 
        self.Checkbutton(fmb3a, "LnPr", self.LnPr_var, self.LnPr, col=1) 
        self.Checkbutton(fmb3a, "TL", self.TL_var, self.TL, col=2) 
        
        tk.Label(fmb3a, text="Invariant-Gamma").grid(row=0, column=3, sticky=tk.W)
        self.Checkbutton(fmb3a, "alpha", self.alpha_var, self.alpha, col=4) 
        self.Checkbutton(fmb3a, "pinvar", self.pinvar_var, self.pinvar, col=5) 
        
        tk.Label(fmb3a, text="Covarion model").grid(row=0, column=6, sticky=tk.W)
        self.Checkbutton(fmb3a, "off-on", self.off_on_var, self.off_on, col=7) 
        self.Checkbutton(fmb3a, "on-off", self.on_off_var, self.on_off, col=8) 
        
        fmb3a.grid(row=row, column=1, sticky=tk.W)
        row += 1

        row = self.break_line( inFrame,  row)

        row = self.Entry(inFrame, 'Change01:', self.mrBayes_replace01_var, self.mrBayes_replace01, row, 100)
        row = self.Entry(inFrame, 'Change02:', self.mrBayes_replace02_var, self.mrBayes_replace02, row, 100)
        row = self.Entry(inFrame, 'Change03:', self.mrBayes_replace03_var, self.mrBayes_replace03, row, 100)
        row = self.Entry(inFrame, 'Change04:', self.mrBayes_replace04_var, self.mrBayes_replace04, row, 100)
        row = self.Entry(inFrame, 'Drop text:', self.string_kill_var, self.string_kill, row, 100)
        row = self.Entry(inFrame, 'Burnin:', self.burnin_var, self.burnin, row, 10)

        row = self.break_line( inFrame,  row)

        f2 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        row = self.ComboEntry(f2, 'Image: ', self.listImage, self.imgType_var, self.imgType, row)
        row = self.Entry2(f2, 'DPI: ', self.dpi_var, self.dpi, row, col=3, width=5)
        f2.grid(row=row, column=1, sticky=tk.W)
        row += 1

        f3 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        self.Checkbutton( f3, "Save data", self.saveData_var, self.saveData, col=0)        
        self.Checkbutton( f3, "See image", self.showGraph_var, self.showGraph, col=1)        
        self.Checkbutton( f3, "Save image", self.saveGraph_var, self.saveGraph, col=2)        
        f3.grid(row=row, column=1, sticky=tk.W)
        row += 1

        row = self.break_line( inFrame,  row)

        row = self.Entry(inFrame, 'Message:', self.msg_obs_var, self.msg_obs, row, lenInputFrame10)


    def form_classifier(self, num, tab):
        # inFrame, row = self.set_inFrame(tab)
        inFrame = tab
        row = 0
        row = self.entry_main_label( inFrame,  row, num)

        row = self.Entry(inFrame, 'Organism:', self.organism_var, self.organism, row,lenInputFrame70,readonly=True)
        row = self.Entry(inFrame, 'Gene:', self.gene_var, self.gene, row,lenInputFrame70,readonly=True)
        row = self.Entry(inFrame, 'Title:', self.title_var, self.title, row,lenInputFrame10,readonly=True)

        #row = self.Choose_MinMax(inFrame, row)
        f1 = tk.Frame(inFrame, width=60, borderwidth=3, relief=tk.RAISED)
        self.seqType_var.set(self.seqType) # initialize
        tk.Radiobutton(f1, text="min", variable=self.minmax_var, value="mincut", indicatoron=0, command = self.on_change_fields).grid(row=0, column=2)
        tk.Radiobutton(f1, text="max", variable=self.minmax_var, value="maxmer", indicatoron=0, command = self.on_change_fields).grid(row=0, column=3)
        f1.grid(row=row, column=0, sticky=tk.W)

        f1b = tk.Frame(inFrame, width=60, borderwidth=3, relief=tk.RAISED)
        self.Checkbutton(f1b, "bias correction", self.withCorrection_var, self.withCorrection, col=0, command=self.on_change_fields)        
        self.Checkbutton(f1b, "mnat", self.mnat_var, self.mnat, col=1)
        f1b.grid(row=row, column=1, sticky=tk.W)
        row += 1

        fr = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        tk.Radiobutton(fr, text="standard", variable=self.rand_method_var, value="standard", indicatoron=0, command = self.on_change_fields).grid(row=0, column=1)
        tk.Radiobutton(fr, text="shuffling", variable=self.rand_method_var, value="shuffle", indicatoron=0, command = self.on_change_fields).grid(row=0, column=2)
        tk.Radiobutton(fr, text="randomize", variable=self.rand_method_var, value="random", indicatoron=0, command = self.on_change_fields).grid(row=0, column=3)
        fr.grid(row=row, column=1, sticky=tk.W)
        row += 1  


        row = self.break_line( inFrame,  row)
        
        f0 = tk.Frame(inFrame, borderwidth=1, relief=tk.RAISED)
        self.EntryFrame(f0, 'Num of letters:', self.numOfLetters_var, self.numOfLetters, row=0, col=0, width=6)
        self.EntryFrame(f0, 'Frame:', self.frame_var, self.frame, row=0, col=2, width=6)
        self.EntryFrame(f0, 'Offset:', self.offset_var, self.offset, row=1, col=0, width=6)
        self.EntryFrame(f0, 'Cutoff:', self.cutoffNumSeq_var, self.cutoffNumSeq, row=1,col=2, width=6)
        self.EntryFrame(f0, 'Min.Samples:', self.NminSamples_var, self.NminSamples, row=2,col=0, width=6)
        self.EntryFrame(f0, 'Num of Loops:', self.numSimLoops_var, self.numSimLoops, row=2,col=2, width=6)
        f0.grid(row=row,column=1, sticky=tk.W)
        row += 1


        row = self.entry_seqType(inFrame, row, False, False) 

        row = self.break_line( inFrame,  row)

        f2 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        self.Checkbutton(f2, "Save data", self.saveData_var, self.saveData, col=0)        
        self.Checkbutton(f2, "See image", self.showGraph_var, self.showGraph, col=1)        
        self.Checkbutton(f2, "Save image", self.saveGraph_var, self.saveGraph, col=2)        
        self.Checkbutton(f2, "Recalculate", self.de_novo_var, self.de_novo, col=3) 
        self.ComboEntry(f2, 'Image: ', self.listImage, self.imgType_var, self.imgType, row=1)
        self.EntryRow(f2, 'DPI: ', self.dpi_var, self.dpi, row=1, col=3, width=5)
        f2.grid(row=row, column=1, sticky=tk.W)
        row += 1

        row = self.break_line( inFrame,  row)
        
        f3 = tk.Frame(inFrame, borderwidth=3, relief=tk.RAISED)
        self.Checkbutton( f3, "Show Message", self.showmessage_var, self.showmessage, col=1)        
        f3.grid(row=row, column=1, sticky=tk.W)
        row += 1

        row = self.Entry(inFrame, 'Message:', self.msg_obs_var, self.msg_obs, row,lenInputFrame10)

        
    def rb_seqType_selected(self):
        self.on_change_fields(None)

    def rb_vert_horiz_selected(self):
        
        if self.vert_horiz_var.get() == 'HMI':
            self.which = "HMI"
            self.input_jsd = self.input_jsd_horiz
        elif self.vert_horiz_var.get() == 'VMI':
            self.which = "VMI"
            self.input_jsd = self.input_jsd_vert
        else:
            self.which = "VSH"
            self.input_jsd = self.input_jsd_vert_entropy
        
        self.unit = 'mnat' if self.mnat else 'nat'
        self.set_cluster_filenames()
        self.input_jsd_var.set(self.input_jsd)
        self.cluster_input_filename_var.set(self.cluster_input_filename)

    
    def rb_seq_yes_no_selected(self):
        # tkMessageBox.showinfo('', str(self.superContig_var.get() ))
        pass


    
    def read_fasta_append_seqs(self, root_filename_fasta, showmessage=False):
        
        try:
            SeqIO.parse(root_filename_fasta, format="fasta").next()
        except:
            stri = 'Could not parse %s'%(root_filename_fasta)
            tkMessageBox.showinfo('Error', stri)
            return False

        for seq_record in SeqIO.parse(root_filename_fasta, format="fasta"):
            ''' maximum merged append '''
            self.sequence.append(SeqRecord(seq_record.seq, id=seq_record.id, description=''))

                    
        return True
    
    def save_fasta_new(self, filename, seqs):
        try:
            handle = open(filename, "w")
            SeqIO.write(seqs, handle, 'fasta')
            ret = True
            stri = 'Saved in %s'%(filename)
            tkMessageBox.showinfo('Warning', stri)
            
        except ValueError:
            stri = "File '%s' not saved. Writing error: %s"%(filename, ValueError)
            tkMessageBox.showinfo('Error', stri)
            ret = False

        finally:   
            handle.close()
            
        return ret
        
    def build_see_tree(self):
        self.get_params()
        self.mb.build_see_tree()

    def build_fig_tree(self):
        self.get_params()
        self.mb.build_fig_tree()

    def run_algor(self):
   
        organism = self.organism_var.get().strip()
        gene = self.gene_var.get().strip()
        title = self.title_var.get().strip()
        
        if organism <> self.organism2 or gene <> self.gene2 or title <> self.title2:
            self.new_search = True

        self.get_params()

        if self.new_search:
            self.saveParams()
        else:
            self.set_paths()
           
        cmp_stri = ''
        

        tk_root.config(cursor="watch")
        
        pipe = None
        
        if self.tab_caption == 'NCBI':
            print 'running pipe01_Find_ManySpecies_OneGene_GBK'
            # thread = Thread(target=pipe1.Pipe(self))
            pipe = pipe1.Pipe(self)
            cmp_stri = 'gbk created: %s\n'%(pipe.filenameGBK)
            cmp_stri += 'logGbkFilename saved %s'%(self.logGbkFilename)
            

        elif self.tab_caption == 'Gbk-Fasta':
            print 'running pipe02_Read_GBK_to_Fasta'
            # thread = Thread(target=pipe2.Pipe(self))
            pipe = pipe2.Pipe(self)
            cmp_stri = 'Fasta files created. Perform alignment to continue the analysis.'
            
            if not pipe.failure:
                self.filter_listSpecies(mode=self.rb_seq_filter_var.get(), cutoff=self.cutoffNumSeq_var.get())
            

        elif self.tab_caption == 'Alignment':
            print 'aligning ..'
            # thread = Thread(target=pipe2.Pipe(self))
            pipe = pipe3.Pipe(self)
            cmp_stri = "Aligned fasta files created, _aligned, mincut and maxmer"
            
            if not pipe.failure:
                self.filter_listSpecies(mode=self.rb_seq_filter_var.get(), cutoff=self.cutoffNumSeq_var.get())
                                        
        elif self.tab_caption == 'Purging':
            print 'running pipe04_purge_repeated_sequences'
            # thread = Thread(pipe3.Pipe(self))
            pipe = pipe4.Pipe(self)
            cmp_stri = ''
            
            if not pipe.failure:
                self.filter_listSpecies(mode=self.rb_seq_filter_var.get(), cutoff=self.cutoffNumSeq_var.get())
            
        elif self.tab_caption == 'Consensus':
            print 'running pipe05_Change_Ambiguous_Nuc'
            # thread = Thread(pipe3.Pipe(self))
            pipe = pipe5.Pipe(self)
            cmp_stri = ''

            if not pipe.failure:
                self.filter_listSpecies(mode=self.rb_seq_filter_var.get(), cutoff=self.cutoffNumSeq_var.get())

        elif self.tab_caption == 'VMI':
            print 'running pipe06_VMI_Distribution_and_JSD'
            pipe = pipe6.Pipe(self, 'Entropy')
            cmp_stri = 'Vertical Mutual Information analyzed.'
            
        elif self.tab_caption == 'HMI':
            print 'running pipe07_Horizontal MI Distribution'
            pipe = pipe7.Pipe(self)
            cmp_stri = 'Horizontal Mutual Information analyzed.'
           
        elif self.tab_caption == 'JSD':
            print 'running pipe09a_JSD'
            pipe = pipe8.Pipe(self)
            cmp_stri = 'JSD analyzed.'
                        
        elif self.tab_caption == 'Cluster':
            print 'running hierarchical cluster'
            pipe = pipe9.Pipe(self)
            cmp_stri = 'end hierarchical cluster.'

        elif self.tab_caption == 'Entropy':
            print 'running Shannon Entropy simulation'
            pipe = pipe10.Pipe(self)
            cmp_stri = 'end Shannon Entropy simulation'

        elif self.tab_caption == 'MrBayes':
            print 'running MrBayes params analyzer'
            pipe = pipe12.Pipe(self)
            cmp_stri = 'end MrBayes params analyze.'

        elif self.tab_caption == 'Classifier':
            print 'running HMI Classifier'
            pipe = pipe13.Pipe(self)
            cmp_stri = 'end HMI Classifier'

        else:
            tkMessageBox.showinfo('Error', 'algorithm not found %s'%(self.tab_caption))
            return
   
        try:
            tk_root.config(cursor='arrow')
        except:
            try:
                tk_root.config(cursor=self.cursor)
            except:
                pass

        if pipe.failure:
            tkMessageBox.showinfo('error', pipe.error_msg)
        else:
            if not self.showGraph and cmp_stri != "":
                tkMessageBox.showinfo('message', 'Task ended.' + cmp_stri)
                
        del pipe
        try:
            gc.collect() 
        except:
            pass        
        
    def run_algor2(self):
        if self.new_search:
            self.saveParams()
                    
        if self.tab_caption == 'VMI':
            self.run_vmi_entropy()
        
    def run_vmi_entropy(self):
        self.set_paths()

        self.cursor = tk_root.cget("cursor")
        tk_root.config(cursor="watch")

        #pipe = None
        print 'running pipe06_Vertical_MI_HeatMap'
        pipe = pipe6.Pipe(self, 'VMI')
        cmp_stri = 'Vertical Mutual Information analyzed.'

        try:
            tk_root.config(cursor=self.cursor)
        except:
            pass

        if pipe.failure:
            tkMessageBox.showinfo('error', pipe.error_msg)
        else:
            if not self.showGraph:
                tkMessageBox.showinfo('message', 'Task ended...' + cmp_stri)

        
    def saveSpeciesListParams(self):
        text = ''
        
        lista = self.dicParams.keys()
        lista = sorted(lista)
        
        for species in lista:
            mat = self.dicParams[species]
            
            line = "%s,%s,%i,%i,%i,%i,%i,%i"%(mat[0],species,mat[1],mat[2],mat[3],mat[4],mat[5],mat[6])

            text += line + '\n'

                    
        self.totList = self.totList_var.get()
        self.totSeqs = self.totSeqs_var.get()
        self.cutoffNumSeq = self.cutoffNumSeq_var.get()
        '''      #Species               #Seqs'''
        text = str(self.totList) + ',' + str(self.totSeqs) + '\n' + \
               'cutoffNumSeq' + ' ' + str(self.cutoffNumSeq) + '\n' + \
                text
        
        self.write(self.filename_species_list, text)

       
                                
    def clear_text_area(self):
        gc.collect()

        ''' if object  self.text_area.delete(0, END) '''
        self.text_area.delete(1.0, tk.END)
        self.showmsg_obs('')
            


  
    def exit(self):
        ''' Run the app '''
        exit()
        
    '''
    def initialize_mia(self):
        path = self.rootTable
        
        try:
            os.stat(path)
            print 'found %s'%(self.rootTable)
        except:
            if not self.consist_dirs():
                self.set_paths()
                if not self.consist_dirs():
                    print("Problems in subdirectories. Could not save the file.")
                    return
            
            try:
                os.stat(path)
                print 'found %s'%(self.rootTable)

            except:
                raise Exception('Could not find %s, exit.'%(self.rootTable))
        
        self.read_default_ini()
    '''
        
    def read_default_ini(self):
        filename = self.root + self.filename_default

        lista = []
        try:
            f = open(filename, 'r')
            
            lines =  f.read().split('\n')
            [organism, gene, title] = lines[0].split(" - ")
            self.first = lines[0]
            lines.pop(0)
            
            lista = [stri.replace("$"," ") for stri in lines if len(stri)> 5]
            self.new_search = False
            
            self.organism2 = organism
            self.gene2 = gene.strip()
            self.title2 = title.strip()

        except:
            lista = ["Drosophila - Adh - "]
            self.first = lista[0]
            organism = "Drosophila"
            gene = "Adh"
            title = ""
            self.new_search = True
            
            self.organism2 = ""
            self.gene2 = ""
            self.title2 = ""
            
            self.aligned_consensus = "consensus"
            self.standard_covarion = "standard"
            self.rand_method = "standard"

       
        self.search_list = lista
        self.search_org_gene_var.set(self.first)
        
        self.get_default_organism_gene(organism, gene, title)
    

    def get_default_organism_gene(self, organism, gene, title):
        self.organism = organism.strip()
        self.gene     = gene.strip()
        self.title    = title.strip()
        
        self.set_gene_title()

        self.organism_var.set(self.organism)
        self.gene_var.set(self.gene)
        self.title_var.set(self.title)

        self.set_paths()
        
        filename = self.rootTable + self.organism + '_' + self.gene_title + '_default.ini'
        
        try:
            if os.path.isfile(filename):
                f = open(filename, 'r')
                
                lines = f.read().split('\t\n')
               
                try:
                    self.dna_prot = lines[0]
                    self.isProtein = (self.dna_prot != 'DNA')
                    
                    self.seqType = lines[1]
                    #self.organism = lines[2]
                    #self.gene =  lines[3]
                    self.species =  ''

                    self.cutoffLength =  int(lines[5])
                    self.cutoffNumSeq =  int(lines[6])
                    self.numOfLetters =  int(lines[7])
                    self.showGraph =  (lines[8] == 'True')
                    self.saveGraph =  (lines[9] == 'True')
                    self.de_novo =  (lines[10] == 'True')
                    
                    self.rootFasta =  self.home + lines[11]
                    self.rootImage =  self.home + lines[12]
                    self.rootTable = self.home + lines[13]

                    self.geneList   = lines[14]
                    self.specieList = lines[15]
                    self.db         = lines[16]
                    
                    self.completeSeq = lines[17]
                    self.completeGen = lines[18]
                    self.chromosome = lines[19]
                    self.shotGun = lines[20]
                    self.contig = lines[21]
                    self.superContig = lines[22]                 
                                                     
                    self.logGbkFilename = lines[23]
                    self.retmax = int(lines[24])
                    
                    self.email = lines[25]
                except:
                    self.dna_prot = 'DNA'
                    self.isProtein = (self.dna_prot != 'DNA')
                    self.seqType = 'Gene'
                    #self.organism = 'Drosophila'
                    #self.gene =  ''
                    self.species =   ''

                    self.cutoffLength =  400
                    self.cutoffNumSeq =  10
                    self.numOfLetters = 1
                    self.showGraph =  True
                    self.saveGraph =  False
                    self.de_novo =  False
                                            
                    self.geneList   =  ''
                    self.specieList =  ''
                    self.db         =  'nucleotide'
                    
                    self.completeSeq = 'no'
                    self.completeGen = 'no'
                    self.chromosome = 'no'
                    self.shotGun = 'no' == 'True'
                    self.contig = 'no'
                    self.superContig = 'no'    
                             
                                                     
                    self.logGbkFilename = 'log_gbk_%s_%s.txt'%(self.organism, self.gene_title)
                    self.retmax = 10000
                    
                    self.email = 'xxx#xxxx.com'
                
                try:
                    self.maxiSeq = lines[26]
                    self.minSeq  = lines[27]
                except:
                    self.maxiSeq = 100000
                    self.minSeq  = 80
                    
                try:
                    self.stop = (lines[28] == 'True')
                    self.showmessage = (lines[29] == 'True')
                    self.listBadTerms = lines[30]
                except:
                    self.stop = False
                    self.showmessage = False
                    self.listBadTerms = 'synthetic' 

                try:
                    self.confirm_gene = (lines[31] == 'True')
                except:
                    self.confirm_gene = False


                try:
                    self.numOfLetters = int(lines[32])
                    self.frame = int(lines[33])
                    self.isLog = (lines[34] == 'True')
                    self.offset = int(lines[35])
                    
                    
                except:
                    self.numOfLetters = 1
                    self.frame = 0
                    self.isLog = False
                    self.offset = 0
                    
                self.each_all = 'each'
                
                try:
                    self.minmax = lines[36]
                    self.partialCDS = lines[37]
                    self.saveData = (lines[38] == 'True')
                    self.title = lines[39]

                except:
                    self.minmax = 'mincut'
                    self.partialCDS = '---'
                    self.saveData  = True
                    self.title = ''
                
                try:
                    self.numOfExperiments = lines[40] 
                    self.sartAt = lines[41] 
                    self.lengthSim = lines[42]
                    self.cluster_method = lines[43]
                    self.vert_horiz = lines[44]
                    

                except:
                    self.numOfExperiments = 10
                    self.sartAt = 5
                    self.lengthSim  = 500
                    self.cluster_method = 'weighted'
                    self.vert_horiz = 'Horizontal'
                    
                try:
                    self.withCorrection = (lines[45]=='True')
                    self.mnat = (lines[46]=='True')
                except:
                    self.withCorrection = False
                    self.mnat = True
                
                try:
                    self.numOfSEs = int(lines[47])
                except:
                    self.numOfSEs = 2
                    
                try:
                    self.norm = (lines[48]=='True')
                except:
                    self.norm = False

                try:
                    self.imgType = (lines[49])
                    self.dpi = int(lines[50])
                except:
                    self.imgType = 'png'
                    self.dpi = 300
                    
                try:
                    self.colorThreshold = float(lines[51])
                except:
                    self.colorThreshold = 0.025
                    

                try:
                    self.alignment = lines[52]
                    self.horizCutoff = int(lines[53])   # deprecated
                    self.minVertCutoff = int(lines[54])
                    self.maxVertCutoff = int(lines[55])
                except:
                    self.alignment = 'Muscle'
                    self.horizCutoff = 50
                    self.minVertCutoff = 10
                    self.maxVertCutoff = 25
                    
                self.realing_each_species = False
                self.realing_all = False
  
                try:
                    self.heatmap_color = lines[56]
                    self.heatmap_ceil_value = float(lines[57])
                    self.heatmap_ceil = (lines[58]=='True')
                except:
                    self.heatmap_color = 'seismic'
                    self.heatmap_ceil_value = 1100
                    self.heatmap_ceil = False
                    
                try:
                    self.is3D = (lines[59]=='True')
                except:
                    self.is3D = False
                    
                try:
                    self.rootTree = self.home + lines[60]
                except:
                    self.rootTree = self.home + '/data/%s/%s/trees/'%(self.organism, self.gene_title)
                                                   
                try:
                    self.rootEntropy = lines[61]
                    if self.rootEntropy == '':
                        self.rootEntropy = self.home + '/data/entropy/'
                except:
                    self.rootEntropy = self.home + '/data/entropy/'
                                  

                try:
                    self.start_stop_codon = (lines[62] == 'True')
                except:
                    self.start_stop_codon = False

                try:
                    self.cds = lines[63]
                    self.mrna = lines[64]
                except:
                    self.cds = 'yes'
                    self.mrna = 'no'
                             
                try:
                    self.mrBayes_path = lines[65]
                    
                    if self.mrBayes_path == "":
                        if self.isWindows:
                            self.mrBayes_path = "c:/mrbayes/"
                        else:
                            self.mrBayes_path = self.home + '/mrbayes/'                    
                    self.piA = (lines[66] == 'True')
                    self.piC = (lines[67] == 'True')
                    self.piG = (lines[68] == 'True')
                    self.piT = (lines[69] == 'True')
                    self.rAC = (lines[70] == 'True')
                    self.rAG = (lines[71] == 'True')
                    self.rAT = (lines[72] == 'True')
                    self.rCG = (lines[73] == 'True')
                    self.rCT = (lines[74] == 'True')
                    self.rGT = (lines[75] == 'True')
                    self.LnL = (lines[76] == 'True')
                    self.LnPr  = (lines[77] == 'True')
                    self.TL    = (lines[78] == 'True')
                    self.alpha = (lines[79] == 'True')
                    self.off_on = (lines[80] == 'True')
                    self.on_off = (lines[81] == 'True')
                    self.pinvar = (lines[82] == 'True')
                    self.LnL_HMean = (lines[83] == 'True')

                except:
                    if self.isWindows:
                        self.mrBayes_path = "c:/mrbayes/"
                    else:
                        self.mrBayes_path = self.home + '/mrbayes/'
        
                    self.piA = True
                    self.piC = True
                    self.piG = True
                    self.piT = True
                    self.rAC = True
                    self.rAG = True
                    self.rAT = True
                    self.rCG = True
                    self.rCT = True
                    self.rGT = True
                    self.LnL = True
                    self.LnPr  = True
                    self.TL    = True
                    self.alpha = True
                    self.off_on = True
                    self.on_off = True
                    self.pinvar = True
                    self.LnL_HMean = True
                                   
                try:
                    self.NminSamples = int(lines[84])
                    self.numSimLoops = int(lines[85])
                except:
                    self.NminSamples = 6
                    self.numSimLoops = 10
                  
                try:
                    self.figtree_path = lines[86]
                    self.mrBayes_replace01 = lines[87]
                    self.mrBayes_replace02 = lines[88]
                    self.mrBayes_replace03 = lines[89]
                    self.mrBayes_replace04 = lines[90]
                    self.string_kill = lines[91]     # '_Adh_Prod_
                    self.colaps_tree = lines[92]=="True"
                    self.aligned_consensus = lines[93]
                except:
                    if self.isWindows:
                        self.figtree_path = "C:/FigTree v1.4.2/FigTree v1.4.2.exe"
                    else:
                        self.figtree_path = "Figtree"
                                
                    self.mrBayes_replace01 = "americana: americana"
                    self.mrBayes_replace02 = "pseudoobscura, persimilis: pseudo_persi"
                    self.mrBayes_replace03 = ""
                    self.mrBayes_replace04 = ""
                    self.string_kill = '_%s_Prod_'%self.gene()
                    self.colaps_tree = True
                    self.aligned_consensus = "consensus"

                try:
                    self.standard_covarion = lines[94]
                    if self.standard_covarion == "":
                        self.standard_covarion = "standard"
                except:
                    self.standard_covarion = "standard"

                try:
                    self.mrBayes_fn01 = lines[95]
                    self.mrBayes_fn02 = lines[96]
                    self.mrBayes_fn03 = lines[97]
                    self.mrBayes_fn04 = lines[98]
                except:
                    ''' mrbayes 
                        Drosophila_maxmer_Gene_Adh_100L_cutoff7_aligned.nxS'''

                    self.mrBayes_fn01 = ('%s_%s_%s_%s_%iL_cutoff%i_%s.nxs') %\
                         (self.organism, "mincut", self.seqType, self.gene_title, self.cutoffLength, self.cutoffNumSeq,"aligned")
                    self.mrBayes_fn02 = ('%s_%s_%s_%s_%iL_cutoff%i_%s.nxs') %\
                         (self.organism, "mincut", self.seqType, self.gene_title, self.cutoffLength, self.cutoffNumSeq,"consensus")
                    self.mrBayes_fn03 = ('%s_%s_%s_%s_%iL_cutoff%i_%s.nxs') %\
                         (self.organism, "maxmer", self.seqType, self.gene_title, self.cutoffLength, self.cutoffNumSeq,"aligned")
                    self.mrBayes_fn04 = ('%s_%s_%s_%s_%iL_cutoff%i_%s.nxs') %\
                         (self.organism, "maxmer", self.seqType, self.gene_title, self.cutoffLength, self.cutoffNumSeq,"consensus")

                try:
                    self.burnin = float(lines[99])
                except:
                    self.burnin = .25
                                       
                try:
                    self.mrBayes_fn01_cov = lines[100]
                    self.mrBayes_fn02_cov = lines[101]
                    self.mrBayes_fn03_cov = lines[102]
                    self.mrBayes_fn04_cov = lines[103]
                except:
                    ''' mrbayes 
                        Drosophila_maxmer_Gene_Adh_100L_cutoff7_aligned.nxS'''

                    self.mrBayes_fn01_cov = ('%s_%s_%s_%s_%iL_cutoff%i_%s_covarion.nxs') %\
                         (self.organism, "mincut", self.seqType, self.gene_title, self.cutoffLength, self.cutoffNumSeq,"aligned")
                    self.mrBayes_fn02_cov = ('%s_%s_%s_%s_%iL_cutoff%i_%s_covarion.nxs') %\
                         (self.organism, "mincut", self.seqType, self.gene_title, self.cutoffLength, self.cutoffNumSeq,"consensus")
                    self.mrBayes_fn03_cov = ('%s_%s_%s_%s_%iL_cutoff%i_%s_covarion.nxs') %\
                         (self.organism, "maxmer", self.seqType, self.gene_title, self.cutoffLength, self.cutoffNumSeq,"aligned")
                    self.mrBayes_fn04_cov = ('%s_%s_%s_%s_%iL_cutoff%i_%s_covarion.nxs') %\
                         (self.organism, "maxmer", self.seqType, self.gene_title, self.cutoffLength, self.cutoffNumSeq,"consensus")

                                                 
                self.rand_method = "standard"            
                self.totList = 0
                
                f.close()
                print '\n-------------------------------'
                print '%s is ok.'%(filename)
                fileError = False
                
            else:
                fileError = True
        except:
            fileError = True

        if fileError and self.first_time:
            print '\n-------------------------------'
            print 'Could not read from %s, default is defined'%(filename)
                
            self.dna_prot = 'DNA'
            self.isProtein = (self.dna_prot != 'DNA')
            self.seqType = 'Gene'
            
            '''  Capsella bursa-pastoris (206), Calystegia soldanella (91), Staphylococcus aureus (83)
                 Betula humilis (44), Arabidopsis thaliana (35), Populus tremula (34), Echinacea angustifolia (24),
                 Echinacea pallida (20), Echinacea laevigata (20), Dioscorea tokoro (19),
                 Pseudomonas aeruginosa (17), Capsella grandiflora (16), Echinacea paradoxa (16),
                 Pyrus pyrifolia (15), Echinacea tennesseensis (14)
            '''
            self.species = ''
            self.cutoffLength = 100
            self.cutoffNumSeq = 10
            self.numOfLetters =  1
            
            self.geneList   =  '' # 'Adh,Adh-1,Adh-2,Adh1,Adh2,AdhPsi,Adhr,Adh-r,AdhT'
            #stri = str(drosGene.speciesList).replace("'","")
            #stri = stri[1:len(stri)-1]
            self.specieList = ''
            self.db = 'Nucleotide'
            
            self.cds = 'yes'
            self.mrna = 'no'
            self.completeGen = 'no'
            self.completeSeq = 'no'
            
            self.chromosome = 'no'
            self.shotGun = 'no'
            self.contig = 'no'
            self.superContig = 'no'
            self.minmax = 'mincut'
            self.each_all = 'each'
            self.partialCDS = '---'
            self.listBadTerms = 'synthetic'
         
            self.retmax = 10000
            self.logGbkFilename = 'log_gbk_%s_%s.txt'%(self.organism, self.gene_title)

            self.email = 'xxx@xxx.com'
            
            self.maxiSeq = 100000
            self.minSeq  = 80
            
            self.stop = False
            self.showmessage = False  
            self.confirm_gene = False    

            self.numOfLetters = 1
            self.frame = 0
            self.isLog = False
            self.showGraph = True
            self.saveGraph = True
            self.offset = 0     
            self.saveData = True
            
            self.numOfExperiments = 10
            self.sartAt = 5
            self.lengthSim  = 500
            
            self.cluster_method = 'weighted'
            self.vert_horiz = 'Horizontal'
            
            self.withCorrection = False
            self.mnat = True
            
            self.de_novo = False
            self.numOfSEs = 2
            self.norm = False
            
            self.imgType = 'png'
            self.dpi = 300
            self.colorThreshold = 0.025

            self.alignment = 'Muscle'
            self.horizCutoff = 50
            self.minVertCutoff = 10
            self.maxVertCutoff = 25
                    
            self.realing_each_species = False
            self.realing_all = False
      
            self.heatmap_color = 'seismic'
            self.heatmap_ceil_value = 1100
            self.heatmap_ceil = False
            self.is3D = False
                    
            self.start_stop_codon = False

        self.set_paths()

        if not self.consist_dirs():
            print("Problems in subdirectories. Could not save the file.")

        self.load_vars()
        
        if self.first_time:
            self.first_time = False
        '''
        else:
            self.build_listSpecies()
            self.filter_listSpecies(mode=self.rb_seq_filter_var.get(), cutoff=self.cutoffNumSeq_var.get())
        '''


        return True
            
    def set_defaults(self, organism='Drosophila', gene='Adh'):
        self.rb_seq_filter = ">="
        self.rb_inv_filter = "inv"
        self.totList = 0
        self.totSeqs = 0
        self.rb_inv_filter = 'inv'
        self.minmax = 'mincut'
        self.each_all = 'each'

        self.organism = organism.strip()
        self.gene = gene.strip()
        self.title = ''

        self.set_gene_title()

        self.organism2 = ""
        self.gene2 = ""
        self.title2 = ""
        
        self.home = os.path.expanduser("~")
        self.home = self.home.replace("\\","/")
        
        self.root = self.home + '/data/'               
        
        self.filename_default = 'default.ini'
        self.shannon_filename = 'shannon_random_DNA_LetterNNN_ExpNNN_dic.txt'
        self.cluster_input_filename = ''
        
        self.mrBayes_replace01 = "americana: americana"
        self.mrBayes_replace02 = "pseudoobscura, persimilis: pseudo_persi"
        self.mrBayes_replace03 = ""
        self.mrBayes_replace04 = ""
        self.string_kill = '_%s_Prod_'%(self.gene_title)
        self.colaps_tree = True
        self.aligned_consensus = "consensus"
        self.standard_covarion = "standard"
        self.rand_method = "standard"

        if self.isWindows:
            self.mrBayes_path = "c:/mrbayes/"
            self.figtree_path = "C:/FigTree v1.4.2/FigTree v1.4.2.exe"
        else:
            self.mrBayes_path = self.home + '/mrbayes/'
            self.figtree_path = "Figtree"
            
        self.which = 'HMI'
        self.stri_random = ""
        self.label_random = ""
        
        self.burnin = .25
                        
                                  
    def set_paths(self):    
        self.set_gene_title()

        self.root_org   = self.home + '/data/%s/' %(self.organism)       
        self.root_files = self.home + '/data/%s/%s/' %(self.organism, self.gene_title)       
        self.rootFasta = self.home + '/data/%s/%s/fasta/'%(self.organism, self.gene_title)
        self.rootImage = self.home + '/data/%s/%s/image/'%(self.organism, self.gene_title)
        self.rootTable = self.home + '/data/%s/%s/files/'%(self.organism, self.gene_title)
        self.rootTree = self.home + '/data/%s/%s/trees/'%(self.organism, self.gene_title)
        self.rootEntropy = self.home + '/data/entropy/'
        
        try:
            self.rootFasta_var.set(self.rootFasta)
            self.rootImage_var.set(self.rootImage)
            self.rootTable_var.set(self.rootTable)
            self.rootTree_var.set(self.rootTree)
            self.rootEntropy_var.set(self.rootEntropy)
        except:
            pass

        if self.isWindows:
            self.root_seaview = "c:/seaview/"
            self.seaview = self.root_seaview + "seaview"
        else:
            self.root_seaview = ""
            self.seaview = "seaview"

    def consist_dirs(self):
        roots = [self.root, self.root_org, self.root_files, self.rootFasta, self.rootImage, self.rootTable, self.rootTree, self.rootEntropy]
        list_roots = ['root', 'root org', 'root Files', 'root Fasta', 'root Image', 'root Table', 'root Trees', 'root Entropy']

        print ''
        for i in range(0, len(roots)):
            root = roots[i]
            ds   = list_roots[i]
            
            path = os.path.dirname(root)         
            try:
                os.stat(path)
                print "%s = %s"%(ds,root)
            except:
                print "creating %s = %s"%(ds,root)
                try:
                    os.mkdir(path)
                except:
                    return False
                
        return True

    def load_vars(self):

        self.organism_var.set(self.organism)
        self.geneList_var.set( self.geneList)
        self.specieList_var.set(self.specieList)
          
        self.species_var.set(self.species)
        
        self.gene = self.gene.strip()
        self.gene_var.set(self.gene)
        
        self.title = self.title.strip()
        self.title_var.set(self.title)
        
        self.seqType_var.set(self.seqType)

        self.dna_prot_var.set(self.dna_prot)
        self.seqType_var.set(self.seqType)

        self.cutoffNumSeq_var.set(self.cutoffNumSeq) 
        self.cutoffLength_var.set(self.cutoffLength)
        self.numOfLetters_var.set(self.numOfLetters)
        
        self.rootFasta_var.set(self.rootFasta)
        self.rootImage_var.set(self.rootImage)
        self.rootTable_var.set(self.rootTable)
        self.rootTree_var.set(self.rootTree)
        self.rootEntropy_var.set(self.rootEntropy)
                                
        self.db_var.set(self.db)
        self.alignment_var.set(self.alignment)
        self.horizCutoff_var.set(self.horizCutoff)
        self.minVertCutoff_var.set(self.minVertCutoff)
        self.maxVertCutoff_var.set(self.maxVertCutoff)

        self.realing_each_species_var.set(False)
        self.realing_all_var.set(False)

        self.cds_var.set(self.cds)
        self.mrna_var.set(self.mrna)
        self.partialCDS_var.set(self.partialCDS)
        self.completeSeq_var.set(self.completeSeq)
        self.completeGen_var.set(self.completeGen)
        self.chromosome_var.set(self.chromosome)
        self.shotGun_var.set(self.shotGun)
        self.contig_var.set(self.contig)
        self.superContig_var.set(self.superContig)
        
        self.minmax_var.set(self.minmax)
        self.each_all_var.set(self.each_all)
        
        self.listBadTerms_var.set(self.listBadTerms)

        self.retmax_var.set(self.retmax)
        self.logGbkFilename_var.set(self.logGbkFilename)
        self.email_var.set(self.email)
        self.maxiSeq_var.set(self.maxiSeq)
        self.minSeq_var.set(self.minSeq)
        
        self.stop_var.set(self.stop)
        self.showmessage_var.set(self.showmessage)
        self.confirm_gene_var.set(self.confirm_gene)
        
        self.numOfLetters_var.set(self.numOfLetters)
        self.frame_var.set(self.frame)
        self.isLog_var.set(self.isLog)
        self.showGraph_var.set(self.showGraph)
        self.saveGraph_var.set(self.saveGraph)
        self.offset_var.set(self.offset)
        self.de_novo_var.set(self.de_novo)

        self.vert_horiz_var.set(self.vert_horiz)
        self.saveData_var.set(self.saveData)

        self.numOfExperiments_var.set(self.numOfExperiments)
        self.sartAt_var.set(self.sartAt)
        self.lengthSim_var.set(self.lengthSim)
        self.shannon_filename_var.set(self.shannon_filename)
        
        
        self.withCorrection_var.set(self.withCorrection)
        self.mnat_var.set(self.mnat)
        self.numOfSEs_var.set(self.numOfSEs)
        self.norm_var.set(self.norm)

        self.cluster_input_filename_var.set(self.cluster_input_filename)
        self.cluster_method_var.set(self.cluster_method)
        self.vert_horiz_var.set(self.vert_horiz_var)

        self.imgType_var.set(self.imgType)
        self.dpi_var.set(self.dpi)
        
        self.colorThreshold_var.set(self.colorThreshold)

        self.heatmap_color_var.set(self.heatmap_color)
        self.heatmap_ceil_value_var.set(self.heatmap_ceil_value)
        self.heatmap_ceil_var.set(self.heatmap_ceil)
        self.is3D_var.set(self.is3D)
                           
        self.set_filenames()
        
        self.aligned_consensus_var.set(self.aligned_consensus)
        self.standard_covarion_var.set(self.standard_covarion)
        self.rand_method_var.set(self.rand_method)


        self.piA_var.set(self.piA)
        self.piC_var.set(self.piC)
        self.piG_var.set(self.piG)
        self.piT_var.set(self.piT)

        self.rAC_var.set(self.rAC)
        self.rAG_var.set(self.rAG)
        self.rAT_var.set(self.rAT)
        self.rCG_var.set(self.rCG)
        self.rCT_var.set(self.rCT)
        self.rGT_var.set(self.rGT)


        self.LnL_var.set(self.LnL)
        self.LnPr_var.set(self.LnPr)
        self.TL_var.set(self.TL)

        self.alpha_var.set(self.alpha)
        self.pinvar_var.set(self.pinvar)
        self.LnL_HMean_var.set(self.LnL_HMean)

        self.off_on_var.set(self.off_on)
        self.on_off_var.set(self.on_off)

        self.NminSamples_var.set(self.NminSamples)
        self.numSimLoops_var.set(self.numSimLoops)

        self.figtree_path_var.set(self.figtree_path)
        self.mrBayes_path_var.set(self.mrBayes_path)
        self.colaps_tree_var.set(self.colaps_tree)
   
        self.mrBayes_replace01_var.set(self.mrBayes_replace01)
        self.mrBayes_replace02_var.set(self.mrBayes_replace02)
        self.mrBayes_replace03_var.set(self.mrBayes_replace03)
        self.mrBayes_replace04_var.set(self.mrBayes_replace04)
        self.string_kill_var.set(self.string_kill)

        self.mrBayes_fn01_var.set(self.mrBayes_fn01)
        self.mrBayes_fn02_var.set(self.mrBayes_fn02)
        self.mrBayes_fn03_var.set(self.mrBayes_fn03)
        self.mrBayes_fn04_var.set(self.mrBayes_fn04)
        self.mrBayes_fn01_cov_var.set(self.mrBayes_fn01_cov)
        self.mrBayes_fn02_cov_var.set(self.mrBayes_fn02_cov)
        self.mrBayes_fn03_cov_var.set(self.mrBayes_fn03_cov)
        self.mrBayes_fn04_cov_var.set(self.mrBayes_fn04_cov)
        
        self.burnin_var.set(self.burnin)
       

    def set_gene_title(self):
        self.gene = self.gene.strip()
        self.title = self.title.strip()

        if self.gene != '':
            self.gene_title = self.gene
            self.graph_gene_title = self.gene
        else:
            self.gene_title = ''
            self.graph_gene_title = ''
        
        if self.title != '':
            if self.gene == '':
                self.gene_title = self.title
                self.graph_gene_title  = self.title
            else:
                self.gene_title += '_' + self.title
                self.graph_gene_title  += ' ' + self.title
                
        if self.gene_title == '':
            self.gene_title = 'undef'
            self.graph_gene_title = 'undef'
        else:
            self.gene_title = self.gene_title.replace(" ","").replace(".","")
            self.graph_gene_title = self.graph_gene_title.replace(" ","").replace(".","")
            
               
    def set_names(self):

        self.organism = self.organism_var.get().strip()
        self.gene = self.gene_var.get().strip()
        self.title = self.title_var.get().strip()   

        self.set_paths()                           
                        
        self.cutoffNumSeq = self.cutoffNumSeq_var.get()
        self.cutoffLength = self.cutoffLength_var.get()
        self.seqType = self.seqType_var.get()
        self.numOfLetters = self.numOfLetters_var.get()
        self.frame = self.frame_var.get()
        
        self.set_filenames()
        
    def set_filenames(self):
        self.set_gene_title() 

        self.withCorrection = self.withCorrection_var.get()
        if self.withCorrection:
            self.filename_correction = '_bias_corr'
        else:
            self.filename_correction = ''

        self.output_filename_gbk = '%s_%s.gbk' % (self.organism, self.gene_title)

        self.output_filename_fasta =  '%s_(species)_(type)_%s_%iL.fasta' % \
             (self.organism, self.gene_title, self.cutoffLength) 

        self.aligned_filename_fasta = '%s_(max/min)_%s_%s_%iL_cutoff%i_aligned.fasta' %\
                 (self.organism, self.seqType, self.gene_title, self.cutoffLength, self.cutoffNumSeq)

        self.mnat = self.mnat_var.get()
        self.unit = 'mnat' if self.mnat else 'nat'                 

        self.set_cluster_filenames()

        self.minmax = self.minmax_var.get()
           
        self.input_filename_align =   ('%s_%s_%s_%s_%iL_cutoff%i_aligned.fasta') % \
             (self.organism, self.minmax, self.seqType, self.gene_title, self.cutoffLength, self.cutoffNumSeq)
                    
        self.output_filename_purged =  ('%s_%s_%s_%s_%iL_cutoff%i_purged.fasta') % \
             (self.organism, self.minmax, self.seqType, self.gene_title, self.cutoffLength, self.cutoffNumSeq)
             
        self.output_filename_consensus = ('%s_%s_%s_%s_%iL_cutoff%i_consensus.fasta') % \
             (self.organism, self.minmax, self.seqType, self.gene_title, self.cutoffLength, self.cutoffNumSeq)

         
        sufix_vert = '%s_%s_(species)_%s_NOL%i_%iL_cutoff%i' %\
        (self.organism, self.minmax, self.seqType, self.numOfLetters, self.cutoffLength, self.cutoffNumSeq)

        sufix_horiz = '%s_%s_(species)_%s_frame%i_NOL%i_%iL_cutoff%i' %\
        (self.organism, self.minmax, self.seqType, self.frame, self.numOfLetters, self.cutoffLength, self.cutoffNumSeq)

        self.output_vmi = 'VMI_%s_mij%s.txt'%(sufix_vert, self.filename_correction)
        self.output_vsh = 'VMI_%s_mij%s.txt'%(sufix_vert, self.filename_correction) 
        self.output_hmi = 'HMI_%s_mij%s.txt'%(sufix_horiz, self.filename_correction)
            
            
        self.input_jsd_vert = self.output_vmi
        self.input_jsd_vert_entropy = self.output_vsh  # vertical shannon entropy
        self.input_jsd_horiz = self.output_hmi
                        
        self.rb_vert_horiz_selected()

        self.logGbkFilename = 'log_gbk_%s_%s.txt'%(self.organism, self.gene_title)

        self.msg_obs = ''
        
        self.filename_species_list =  self.rootTable + self.organism + '_' + self.gene_title + "_species_list.txt"  



    def set_cluster_filenames(self):
        self.set_gene_title()

        if self.withCorrection:
            self.filename_correction = '_bias_corr'
            self.str_correction = '-bias corr.'
        else:
            self.filename_correction = ''
            self.str_correction = ''
            

        if self.frame == 0:
            self.str_frame1 = ''
            self.str_frame2 = ''
        else:
            self.str_frame1 = 'frame=%i_'%(self.frame)
            self.str_frame2 = ' frame=%i, '%(self.frame)

        self.sufixTitle = '%s %s%s %s, gene=%s%s NOL=%i' %\
                (self.organism, self.minmax, self.str_correction, self.seqType, self.graph_gene_title, self.str_frame2, self.numOfLetters)
        self.sufix =  ('%s_%s_%s_%s_%sNOL%i_%iL_cutoff%i%s') %\
              (self.organism, self.minmax, self.seqType, self.gene_title,self.str_frame1,  \
               self.numOfLetters, self.cutoffLength, self.cutoffNumSeq, self.stri_random)

        self.mnat = self.mnat_var.get()
        self.unit = 'mnat' if self.mnat else 'nat'                 

        if self.which == 'HMI':
            self.prefix = 'JSD_HMI_'
            self.title_jsd = 'JSD(Horizontal Mutual Information) %s %s\n%s'%\
                (self.label_random, self.unit, self.sufixTitle)

            self.cluster_input_filename = 'JSD_HMI_%s%s%s.txt' % (self.sufix, self.filename_correction, self.stri_random)
            self.striMatrix = "Horizontal MI distance matrix"
            if self.label_random != "":
                self.striMatrix += " - " + self.label_random
                
        elif self.which == 'VMI':
            self.prefix = 'JSD_VMI_'
            self.title_jsd = 'JSD(Vertical Mutual Information) %s %s\n%s'%\
                (self.label_random, self.unit, self.sufixTitle)

            self.cluster_input_filename = 'JSD_VMI_%s%s%s.txt' % (self.sufix, self.filename_correction, self.stri_random)
            self.striMatrix = "Vertical MI distance matrix"
        else:
            self.prefix = 'JSD_VSH_'
            self.title_jsd = 'JSD(Vertical Shannon Entropy) %s %s\n%s'%\
                (self.label_random, self.unit, self.sufixTitle)

            self.cluster_input_filename = 'JSD_VSH_%s%s%s.txt' % (self.sufix, self.filename_correction, self.stri_random)
            self.striMatrix = "Vertical Shannon Entropy"

        self.cluster_input_filename_var.set(self.cluster_input_filename)
        
    def drawSpeciesList(self):
        self.frameSpeciesList = tk.Frame(borderwidth=3, relief=tk.RAISED)
        self.build_lateral_list(self.cmd_filter_speciesList, self.cmd_inv_speciesList)
        
        self.frameSpeciesList.grid(rowspan=3, row=0, column=1, sticky='NSWE')

    def build_lateral_list(self, command1, command2):
        
        supFrame = tk.Frame(self.frameSpeciesList)
        self.Entry(supFrame,'#Species', self.totList_var, self.totList, 0,5)
        self.Entry(supFrame,'#Seqs', self.totSeqs_var, self.totSeqs, 1,5)
        self.Entry(supFrame ,'Cutoff #Seqs', self.cutoffNumSeq_var, self.cutoffNumSeq, 2,5)
        self.Entry(supFrame, 'Cutoff Length:', self.cutoffLength_var, self.cutoffLength, 3,10)
        supFrame.grid(row=0, column=0, sticky=tk.W)

        f0 = tk.Frame(self.frameSpeciesList)
        f1 = tk.Frame(f0, borderwidth=3, width=50, relief=tk.RAISED)
        self.rb_seq_filter_var.set(self.rb_seq_filter)
    
        tk.Radiobutton(f1, text=">=", variable=self.rb_seq_filter_var, value=">=", indicatoron=0, command=command1).grid(row=0, column=0)
        tk.Radiobutton(f1, text="all", variable=self.rb_seq_filter_var, value="all", indicatoron=0, command=command1).grid(row=0, column=1)
        tk.Radiobutton(f1, text="<", variable=self.rb_seq_filter_var, value="<", indicatoron=0, command=command1).grid(row=0, column=2)
        
        f2 = tk.Frame(f0, borderwidth=3, relief=tk.RAISED)
        self.rb_inv_filter_var.set(self.rb_inv_filter)
    
        # tk.Button(frame, text=label, command=command).grid(row=row, column=col, sticky=tk.W)

        tk.Button(f2, text="inv", command=self.invert).grid(row=0, column=0)
        tk.Button(f2, text="clear", command=self.reset_all ).grid(row=0, column=1)
        tk.Button(f2, text="set", command=self.set_all).grid(row=0, column=2)
        tk.Button(f2, text="del", command=self.delete_record ).grid(row=0, column=3)
        
        '''
        tk.Button(f2, text="inv", variable=self.rb_inv_filter_var, value="inv", indicatoron=0, command=command2).grid(row=0, column=0)
        tk.Button(f2, text="clear", variable=self.rb_inv_filter_var, value="clear", indicatoron=0, command=command2).grid(row=0, column=1)
        tk.Button(f2, text="set", variable=self.rb_inv_filter_var, value="set", indicatoron=0, command=command2).grid(row=0, column=2)
        tk.Button(f2, text="del", variable=self.rb_inv_filter_var, value="del", indicatoron=0, command=command2).grid(row=0, column=3)
        '''
        f0.grid(column=0, sticky=tk.W)
                
        f1.grid(row=0, column=0, sticky=tk.W)
        tk.Label(f0, text="      ").grid(row=0, column=1, sticky=tk.W)
        f2.grid(row=0, column=2, sticky=tk.W)


        stri = "     num.seqs.         length        species"
        tk.Label(self.frameSpeciesList, text=stri,font=self.cour10).grid(column=0, sticky=tk.W)
        stri = "    ori/min/max     ori/  min/  max  "
        tk.Label(self.frameSpeciesList, text=stri,font=self.cour10).grid(column=0, sticky=tk.W)
            
        '''--------------- species list -----------------'''
        intFrame = tk.Frame(self.frameSpeciesList)
        self.speciesListbox = tk.Listbox(intFrame, width=self.speciesListWidth, height=self.speciesListHeight, font=self.cour10) # width=speciesListWidth,
        self.speciesListbox.bind('<<ListboxSelect>>', self.onClickSpeciesBox)
        self.speciesListbox.bind("<Double-Button-1>", self.onClickSpeciesBoxInvert)
        
                
        ''' orient="vertical", '''
        scrl = tk.Scrollbar(intFrame, orient="vertical", command=self.speciesListbox.yview) 
        self.speciesListbox.config(yscrollcommand=scrl.set)


        self.speciesListbox.grid(column=0, sticky='NSEW')
        scrl.grid(row=0, column=1, sticky='NSE')

        intFrame.grid(column=0, sticky='NSEW')


    def cmd_filter_speciesList(self):
        self.filter_listSpecies(mode=self.rb_seq_filter_var.get(), cutoff=self.cutoffNumSeq_var.get())


    def cmd_search_list(self, event):
        
        self.first = self.search_org_gene_var.get()
        [organism, gene, title] = self.first.split(" - ")
        self.get_default_organism_gene(organism, gene, title)

        self.new_search = False

        self.build_listSpecies()
        self.filter_listSpecies(mode=self.rb_seq_filter_var.get(), cutoff=self.cutoffNumSeq_var.get())

    def invert(self):
        for spec in self.dicParams.keys():
            mat = self.dicParams[spec]
            if mat[0] == 'x':
                mat[0] = ' '
            else:
                mat[0] = 'x'
                
        self.filter_listSpecies(mode=self.rb_seq_filter_var.get(), cutoff=self.cutoffNumSeq_var.get())

        
    def delete_record(self):
        '''
        lista = self.speciesListbox
        index = int(lista.curselection()[0])
        self.delete_record(index)
        '''
        if not tkMessageBox.askyesno("Delete", "Do you want to delete this record?"):
            return
                        
        index = int(self.speciesListbox.curselection()[0])

        species = self.find_species(self.speciesListbox.get(index))
        self.dicParams.pop(species)
        
        self.filter_listSpecies(mode=self.rb_seq_filter_var.get(), cutoff=self.cutoffNumSeq_var.get())
        
        tkMessageBox.showinfo("Warning", "At the end you must save to effectively delete the record.")

    def reset_all(self):
        for spec in self.dicParams.keys():
            self.dicParams[spec][0] = ' '
                
        self.filter_listSpecies(mode=self.rb_seq_filter_var.get(), cutoff=self.cutoffNumSeq_var.get())

    def set_all(self):
        for spec in self.dicParams.keys():
            self.dicParams[spec][0] = 'x'
                
        self.filter_listSpecies(mode=self.rb_seq_filter_var.get(), cutoff=self.cutoffNumSeq_var.get())

        

    def onClickSpeciesBox(self, e):
        stri = self.speciesListbox.selection_get()
        self.species = self.find_species(stri)
        self.species_var.set(self.species)

        
    def onClickSpeciesBoxInvert(self, evt):
        lista = self.speciesListbox
        index = int(lista.curselection()[0])
        stri = lista.get(index)
        
        species = self.find_species(stri)

        check = stri[:3]
        stri = stri[3:]
        
        # print self.speciesListbox.selection_get()
        if check == '[x]':
            stri = '[ ]' + stri
            self.dicParams[species][0] = ' '
        else:
            stri = '[x]' + stri
            self.dicParams[species][0] = 'x'
            
        self.speciesListbox.delete(index, index)
        self.speciesListbox.insert(index, stri)
        self.speciesListbox.index(index)
            


        
    def cmd_inv_speciesList(self):
        if self.rb_inv_filter_var.get() == "inv":
            self.invert()
        elif self.rb_inv_filter_var.get() == "del":
            self.delete_record()
        elif self.rb_inv_filter_var.get() == "clear":
            self.reset_all()
        else:
            self.set_all()
            
    def build_listSpecies(self):
       
        lista = self.read_text_to_species_list(self.filename_species_list)
        self.dicParams = {}
        
        
        if not lista:
            self.totList = 0
            self.totSeqs = 0

            
            try:
                self.speciesListbox.delete(0, tk.END)
            except:
                pass
            return

        ''' two initial dummy lines with totals '''
        try:
            
            stri = lista[0]
            stri = stri.split(',')
            
            stri = lista[1]
            stri = stri.split(' ')
            self.cutoffNumSeq = int(stri[1])
            
        except:
            self.totList = 0
            self.totSeqs = 0
    
        finally:
            try:
                lista.pop(0)
                lista.pop(0)
            except:
                pass
    
        for line in lista:
            mat2 = line.split(",")
            species = mat2[1]
            ''' mat[0] = 'x' '''
            for i in range(2,8):
                mat2[i] = int(mat2[i])
                
            self.dicParams[species] = [mat2[0], mat2[2],mat2[3],mat2[4],mat2[5],mat2[6],mat2[7]]
            # mat = self.dicParams[species]
        


    def filter_listSpecies(self, mode="all", cutoff=10, showmessage=False):
        
        if not self.dicParams:
            self.totList = 0
            self.totSeqs = 0
            return

        self.speciesList = []
        self.totList = 0
        self.totSeqs = 0

        specList = self.dicParams.keys()
        specList = sorted(specList)
                    
        if mode=="all":
            for species in specList:
                mat = self.dicParams[species]
                self.speciesList.append(self.mask%(mat[0],mat[1],mat[2],mat[3],mat[4],mat[5],mat[6], species))   
                self.totList += 1
                self.totSeqs += mat[2]  
        else:            
            if mode==">=":
                for species in specList:
                    mat = self.dicParams[species]
                    if mat[1] >= cutoff:
                        self.speciesList.append(self.mask%(mat[0],mat[1],mat[2],mat[3],mat[4],mat[5],mat[6], species))   
                        self.totList += 1
                        self.totSeqs += mat[2]                     
            else:         
                for species in specList:
                    mat = self.dicParams[species]
                    if mat[1] < cutoff:
                        self.speciesList.append(self.mask%(mat[0],mat[1],mat[2],mat[3],mat[4],mat[5],mat[6], species))   
                        self.totList += 1
                        self.totSeqs += mat[2]  
        

        self.totList_var.set(self.totList)
        self.totSeqs_var.set(self.totSeqs)

        self.speciesListbox.delete(0, tk.END)

        for stri in self.speciesList:
            self.speciesListbox.insert(tk.END, stri)
            
            if showmessage:
                self.showmsg_obs(stri)
                        

    def drawText(self, row):
        frameFooter = tk.Frame(self.tk_root, borderwidth=3, relief=tk.RAISED)

        ''' must review !!! '''       
        if self.isWindows:
            textAreaHeight = 10
            textAreaWidth = 90
        else:
            textAreaHeight = 30
            textAreaWidth = 148

        
        self.text_area = tk.Text(frameFooter,  height=textAreaHeight, width=textAreaWidth,  bg='light cyan') # 
 
        ''' orient="vertical", '''
        scrl = tk.Scrollbar(frameFooter, orient="vertical", command=self.text_area.yview)
        self.text_area.config(yscrollcommand=scrl.set)
               
        sys.stdout = StdoutDirector(self.text_area)

        self.text_area.grid(row=0, column=0, sticky='NSW') # row=row,column=0,
        scrl.grid(row=0, column=1, sticky='NSE') # row=0, orient='vertical', 
        
        frameFooter.grid(columnspan=3, row=row, column=0, sticky='NSWE')  # columnspan=2,  row=row, column=0

        self.text_area.columnconfigure(0, weight=1)

        ''''
        self.nb.update()
        print self.nb.winfo_width() 
        print frameFooter.winfo_width() 
        print self.text_area.winfo_width()
        '''

    
          
    def drawButtons(self, row):
        self.framBtn = tk.Frame(borderwidth=3, relief=tk.RAISED)        
        
        ''' The tk.Button without grid implicit, otherwise, will not be referenced '''
        self.buttonRun1  = tk.Button(self.framBtn, text='Run',    width=24, command=self.run_algor)
        self.buttonRun1.grid(row=0, column=0)
        self.buttonRun2  = tk.Button(self.framBtn, text='Run2',   width=24, command=self.run_algor2)
        self.buttonRun2.grid(row=0, column=1)
        self.buttonSave  = tk.Button(self.framBtn, text ="Save",  width=15, command=self.saveParams)
        self.buttonSave.grid(row=0,column=2)
        self.buttonExit  = tk.Button(self.framBtn, text ="Exit",  width=15, command=self.exit)
        self.buttonExit.grid(row=0,column=3)
        self.buttonClear = tk.Button(self.framBtn, text ="Clear", width=15, command=self.clear_text_area)
        self.buttonClear.grid(row=0,column=4)

        self.framBtn.grid(row=row, sticky='SWE')

    def tabChangedEvent(self,event):
        # print event.widget.tab(event.widget.index("current"),"text")
        #self.labelroot_var.set(event.widget.tab(event.widget.index("current"),"text"))
       
        self.tab_caption = event.widget.tab(event.widget.index("current"),"text")
        
        if self.tab_caption == 'VMI':
            self.buttonRun1['text'] = 'Vertical Entropy'
            self.buttonRun2['text'] = 'Vertical Mutual Information'

            #self.buttonRun2  = tk.Button(self.drawBtn_Frame, text='Vertical Mutual Information',   width=20, command=self.run_algor2)
            #self.buttonRun2.grid(row=0, column=1)
        else:
            # self.buttonRun2.grid_remove()
            try:
                self.buttonRun2['text'] = ''
            except:
                return
            
        if self.tab_caption == 'NCBI':
            self.buttonRun1['text'] = 'Get GBK from NCBI'

        elif self.tab_caption == 'Gbk-Fasta':
            self.buttonRun1['text'] = 'Split Gbk into Fasta'

        elif self.tab_caption == 'Alignment':
            self.buttonRun1['text'] = 'Align'
                           
        elif self.tab_caption == 'Purging':
            self.buttonRun1['text'] = 'Purge'

        elif self.tab_caption == 'Consensus':
            self.buttonRun1['text'] = 'Consensus / Split species'
            
        elif self.tab_caption == 'HMI':
            self.buttonRun1['text'] = 'Horizontal Mutual Information'
           
        elif self.tab_caption == 'JSD':
            self.buttonRun1['text'] = 'Jensen-Shannon Divergence'
                        
        elif self.tab_caption == 'Cluster':
            self.buttonRun1['text'] = 'Hierarchical Cluster'

        elif self.tab_caption == 'Entropy':
            self.buttonRun1['text'] = 'Shannon Entropy'

        elif self.tab_caption == 'MrBayes':
            self.buttonRun1['text'] = 'run params'

        elif self.tab_caption == 'Classifier':
            self.buttonRun1['text'] = 'run Classification'
        
        
    def get_params(self):
        self.organism = self.organism_var.get()

        if self.organism == "Drosophila":
            self.dr = dro.Drosophila()
        else:
            self.dr = None
                    
        self.gene = self.gene_var.get().strip()
        self.geneList = self.geneList_var.get()
        self.specieList = self.specieList_var.get()
        self.species = self.species_var.get()
        self.seqType = self.seqType_var.get()
        self.title = self.title_var.get().strip()
        
        self.set_gene_title()

        self.dna_prot = self.dna_prot_var.get()
        self.isProtein = (self.dna_prot == 'PROT')
        
        self.cutoffNumSeq = self.cutoffNumSeq_var.get()
        self.cutoffLength = self.cutoffLength_var.get()
        self.numOfLetters = self.numOfLetters_var.get()

        self.rootFasta = self.rootFasta_var.get()
        self.rootImage = self.rootImage_var.get()
        self.rootTable = self.rootTable_var.get()
        self.rootTree = self.rootTree_var.get()
        self.rootEntropy = self.rootEntropy_var.get()
    

        self.db = self.db_var.get()
        self.alignment = self.alignment_var.get()
        self.horizCutoff = self.horizCutoff_var.get()
        self.minVertCutoff = self.minVertCutoff_var.get()
        self.maxVertCutoff = self.maxVertCutoff_var.get()
                              

        self.cds = self.cds_var.get()
        self.mrna = self.mrna_var.get()
        self.partialCDS = self.partialCDS_var.get()
        self.completeSeq = self.completeSeq_var.get()
        self.completeGen = self.completeGen_var.get()
        self.chromosome = self.chromosome_var.get()
        self.shotGun = self.shotGun_var.get()
        self.contig = self.contig_var.get()
        self.superContig = self.superContig_var.get()

        self.maxVertCutoff = self.maxVertCutoff_var.get()
        self.minmax = self.minmax_var.get()
        self.each_all = self.each_all_var.get()


        self.listBadTerms = self.listBadTerms_var.get()

        self.retmax = self.retmax_var.get()
        self.logGbkFilename = self.logGbkFilename_var.get()
        self.email = self.email_var.get()
        self.maxiSeq = self.maxiSeq_var.get()
        self.minSeq = self.minSeq_var.get()

        self.stop = self.stop_var.get()
        self.start_stop_codon = self.start_stop_codon_var.get()
        self.showmessage = self.showmessage_var.get()
        self.confirm_gene = self.confirm_gene_var.get()


        self.numOfLetters = self.numOfLetters_var.get()
        self.frame = self.frame_var.get()
        self.isLog = self.isLog_var.get()
        self.showGraph = self.showGraph_var.get()
        self.saveGraph = self.saveGraph_var.get()
        self.saveData = self.saveData_var.get()
        self.offset = self.offset_var.get()
        self.de_novo = self.de_novo_var.get()
        self.imgType = self.imgType_var.get()
        self.dpi = self.dpi_var.get()
       

        self.vert_horiz = self.vert_horiz_var.get()

        if self.vert_horiz ==  'VMI':
            self.which = 'VMI'
        elif self.vert_horiz ==  'VSH':
            self.which = 'VSH'
        else:
            self.which = 'HMI'
                
        self.numOfExperiments = self.numOfExperiments_var.get()
        self.sartAt = self.sartAt_var.get()
        self.lengthSim = self.lengthSim_var.get()
        self.shannon_filename = self.lengthSim_var.get()

        
        self.withCorrection = self.withCorrection_var.get()
        
        if self.withCorrection:
            self.filename_correction = '_bias_corr'
            self.str_correction = '-bias corr.'
        else:
            self.filename_correction = ''
            self.str_correction = ''
                        
        self.mnat = self.mnat_var.get()
        self.numOfSEs = self.numOfSEs_var.get()
        self.norm = self.norm_var.get()

        self.unit = 'mnat' if self.mnat else 'nat'
        self.factor = 1000 if self.mnat else 1
        
        self.cluster_input_filename = self.cluster_input_filename_var.get()
        self.cluster_method = self.cluster_method_var.get()

 
        self.colorThreshold = self.colorThreshold_var.get()
        self.heatmap_color = self.heatmap_color_var.get()
        self.heatmap_ceil_value = self.heatmap_ceil_value_var.get()
        self.heatmap_ceil = self.heatmap_ceil_var.get()
        self.is3D = self.is3D_var.get()
                          
        self.aligned_consensus = self.aligned_consensus_var.get()
        self.standard_covarion = self.standard_covarion_var.get()
        self.rand_method = self.rand_method_var.get()

        self.piA = self.piA_var.get()
        self.piC = self.piC_var.get()
        self.piG = self.piG_var.get()
        self.piT = self.piT_var.get()


        self.rAC = self.rAC_var.get()
        self.rAG = self.rAG_var.get()
        self.rAT = self.rAT_var.get()
        self.rCG = self.rCG_var.get()
        self.rCT = self.rCT_var.get()
        self.rGT = self.rGT_var.get()


        self.LnL = self.LnL_var.get()
        self.LnPr = self.LnPr_var.get()
        self.TL = self.TL_var.get()
        self.alpha = self.alpha_var.get()
        self.pinvar = self.pinvar_var.get()
        self.LnL_HMean = self.LnL_HMean_var.get()


        self.off_on = self.off_on_var.get()
        self.on_off = self.on_off_var.get()

        self.mrBayes_replace01 = self.mrBayes_replace01_var.get()
        self.mrBayes_replace02 = self.mrBayes_replace02_var.get()
        self.mrBayes_replace03 = self.mrBayes_replace03_var.get()
        self.mrBayes_replace04 = self.mrBayes_replace04_var.get()
        self.string_kill = self.string_kill_var.get()

        self.rootFasta = self.rootFasta_var.get()
        self.rootImage = self.rootImage_var.get()
        self.rootTable = self.rootTable_var.get()
        self.rootTree = self.rootTree_var.get()
        self.rootEntropy = self.rootEntropy_var.get()
        self.mrBayes_path = self.mrBayes_path_var.get()
        self.figtree_path = self.figtree_path_var.get()

        self.mrBayes_fn01 = self.mrBayes_fn01_var.get()
        self.mrBayes_fn02 = self.mrBayes_fn02_var.get()
        self.mrBayes_fn03 = self.mrBayes_fn03_var.get()
        self.mrBayes_fn04 = self.mrBayes_fn04_var.get()
        
        self.mrBayes_fn01_cov = self.mrBayes_fn01_cov_var.get()
        self.mrBayes_fn02_cov = self.mrBayes_fn02_cov_var.get()
        self.mrBayes_fn03_cov = self.mrBayes_fn03_cov_var.get()
        self.mrBayes_fn04_cov = self.mrBayes_fn04_cov_var.get()

        self.burnin = self.burnin_var.get()

        self.set_filenames()


    def save_search_list(self):
        self.first = self.organism + " - " + self.gene + " - " + self.title

        lista = copy.deepcopy(self.search_list)
        if self.new_search or self.first not in lista:
            lista.append(self.first)
            
        lista = set(lista)
        lista = list(lista)
        lista = sorted(lista)
        self.search_org_gene_var.set(self.first)        

  
        stri = self.first + '\n'
        for org_gene in lista:
            stri += org_gene + '\n'

        ''' saveing default.ini '''
        filename = self.root + self.filename_default
        self.write(filename, stri)
        
        self.search_list = lista
        self.updtcblist()

                
    def delete_org_gene(self):
        lista = copy.deepcopy(self.search_list)
        lista.remove(self.first)
        
        self.search_org_gene_var.set(lista[0])
        self.first = lista[0]

        stri = self.first + '\n'
        for org_gene in lista:
            stri += org_gene + '\n'

        ''' saveing default.ini '''
        filename = self.root + self.filename_default
        if not self.write(filename, stri):
            return
        
        self.read_default_ini()
        self.updtcblist()
        
        
    def save_all_param_files(self):
        ''' get organism, species ... '''
        self.get_params()
        ''' write self.dicParams. in self.dicParams.'''
        self.saveSpeciesListParams()
        self.save_search_list()
        
                
    def saveParams(self):
        self.set_paths()
        if not self.consist_dirs():
            print("Problems in subdirectories. Could not save the file.")
            return

        self.save_all_param_files()
        self.updtcblist()     

        self.organism2 = self.organism
        self.gene2 = self.gene
        self.title2 = self.title 
        
        if self.db_var.get() == 'Protein':
            self.dna_prot_var.set("PROT")
        else:
            self.dna_prot_var.set("DNA")
            
        self.dna_prot = self.dna_prot_var.get()
        self.isProtein = (self.dna_prot != 'DNA')
            
        stri = ''
        stri += self.dna_prot_var.get() + '\t\n'
        stri += self.seqType_var.get() + '\t\n'
        stri += self.organism_var.get() + '\t\n'
        stri += self.gene_var.get() + '\t\n'
        stri += self.species_var.get() + '\t\n'
        stri += str(self.cutoffLength_var.get()) + '\t\n'
        stri += str(self.cutoffNumSeq_var.get()) + '\t\n'
        stri += str(self.numOfLetters_var.get()) + '\t\n'
        stri += str(self.showGraph_var.get()==True) + '\t\n'
        stri += str(self.saveGraph_var.get()==True) + '\t\n'
        stri += str(self.de_novo_var.get()==True) + '\t\n'

        stri += self.rootFasta_var.get() + '\t\n'
        stri += self.rootImage_var.get() + '\t\n'
        stri += self.rootTable_var.get() + '\t\n'


        stri += self.geneList_var.get() + '\t\n'
        stri += self.specieList_var.get() + '\t\n'
        stri += self.db_var.get()  + '\t\n'
        
        stri += self.completeSeq_var.get()  + '\t\n'
        stri += self.completeGen_var.get()  + '\t\n'
        stri += self.chromosome_var.get()  + '\t\n'
        stri += self.shotGun_var.get()  + '\t\n'
        stri += self.contig_var.get()  + '\t\n'
        stri += self.superContig_var.get()  + '\t\n'
        
        stri += self.logGbkFilename_var.get()  + '\t\n'
        stri += str(self.retmax_var.get())  + '\t\n'
        stri += self.email_var.get() + '\t\n'
        
        stri += str(self.maxiSeq_var.get()) + '\t\n'
        stri += str(self.minSeq_var.get()) + '\t\n'
        stri += str(self.stop_var.get()==True) + '\t\n'
        stri += str(self.showmessage_var.get()==True) + '\t\n'
        stri += self.listBadTerms_var.get() + '\t\n'
        stri += str(self.confirm_gene_var.get()==True) + '\t\n'
        
        stri += str(self.numOfLetters_var.get())  + '\t\n'
        stri += str(self.frame_var.get())  + '\t\n'
        stri += str(self.isLog_var.get()==True) + '\t\n'
        stri += str(self.offset_var.get())  + '\t\n'

        stri += self.minmax_var.get() + '\t\n'
        stri += self.partialCDS_var.get() + '\t\n'
        stri += str(self.saveData_var.get()==True) + '\t\n'
        stri += self.title_var.get() + '\t\n'

        stri += str(self.numOfExperiments_var.get())  + '\t\n'
        stri += str(self.sartAt_var.get())  + '\t\n'
        stri += str(self.lengthSim_var.get())  + '\t\n'

        stri += str(self.cluster_method_var.get())  + '\t\n'
        stri += str(self.vert_horiz_var.get())  + '\t\n'
        
        stri += str(self.withCorrection_var.get()==True)  + '\t\n'
        stri += str(self.mnat_var.get()==True)  + '\t\n'
        stri += str(self.numOfSEs_var.get())  + '\t\n'
        stri += str(self.norm_var.get()==True) + '\t\n'
        stri += self.imgType_var.get() + '\t\n'
        stri += str(self.dpi_var.get()) + '\t\n'
        stri += str(self.colorThreshold_var.get()) + '\t\n'
        
        stri += self.alignment_var.get()  + '\t\n'
        stri += str(self.horizCutoff_var.get())  + '\t\n'
        stri += str(self.minVertCutoff_var.get())  + '\t\n'
        stri += str(self.maxVertCutoff_var.get())  + '\t\n'

        stri += self.heatmap_color_var.get()  + '\t\n'
        stri += str(self.heatmap_ceil_value_var.get())  + '\t\n'
        stri += str(self.heatmap_ceil_var.get()==True)  + '\t\n'
        try:
            stri += str(self.is3D_var.get()==True)  + '\t\n'
        except:
            stri += 'False\t\n'
       
        stri += self.rootTree_var.get() + '\t\n'
        stri += self.rootEntropy_var.get() + '\t\n'
        stri += str(self.start_stop_codon_var.get()==True)  + '\t\n'

        stri += self.cds_var.get()  + '\t\n'
        stri += self.mrna_var.get()  + '\t\n'

        stri += self.mrBayes_path_var.get()  + '\t\n'
        stri += str(self.piA_var.get()==True)  + '\t\n'
        stri += str(self.piC_var.get()==True)  + '\t\n'
        stri += str(self.piG_var.get()==True)  + '\t\n'
        stri += str(self.piT_var.get()==True)  + '\t\n'

        stri += str(self.rAC_var.get()==True)  + '\t\n'
        stri += str(self.rAG_var.get()==True)  + '\t\n'
        stri += str(self.rAT_var.get()==True)  + '\t\n'
        stri += str(self.rCG_var.get()==True)  + '\t\n'
        stri += str(self.rCT_var.get()==True)  + '\t\n'
        stri += str(self.rGT_var.get()==True)  + '\t\n'
        stri += str(self.LnL_var.get()==True)  + '\t\n'
        stri += str(self.LnPr_var.get()==True)  + '\t\n'
        stri += str(self.TL_var.get()==True)  + '\t\n'
        stri += str(self.alpha_var.get()==True)  + '\t\n'
        stri += str(self.off_on_var.get()==True)  + '\t\n'
        stri += str(self.on_off_var.get()==True)  + '\t\n'
        stri += str(self.pinvar_var.get()==True)  + '\t\n'
        stri += str(self.LnL_HMean_var.get()==True)  + '\t\n'
                    
        stri += str(self.NminSamples_var.get())  + '\t\n'
        stri += str(self.numSimLoops_var.get())  + '\t\n'
        
        stri += str(self.figtree_path_var.get())  + '\t\n'
        stri += str(self.mrBayes_replace01_var.get())  + '\t\n'
        stri += str(self.mrBayes_replace02_var.get())  + '\t\n'
        stri += str(self.mrBayes_replace03_var.get())  + '\t\n'
        stri += str(self.mrBayes_replace04_var.get())  + '\t\n' 
        stri += str(self.string_kill_var.get())  + '\t\n'
        stri += str(self.colaps_tree_var.get())  + '\t\n'
        stri += str(self.aligned_consensus_var.get())  + '\t\n'
        stri += str(self.standard_covarion_var.get())  + '\t\n'
  
        stri += str(self.mrBayes_fn01_var.get())  + '\t\n'
        stri += str(self.mrBayes_fn02_var.get())  + '\t\n'
        stri += str(self.mrBayes_fn03_var.get())  + '\t\n'
        stri += str(self.mrBayes_fn04_var.get())  + '\t\n'
     
        stri += str(self.burnin_var.get())  + '\t\n'
        
        stri += str(self.mrBayes_fn01_cov_var.get())  + '\t\n'
        stri += str(self.mrBayes_fn02_cov_var.get())  + '\t\n'
        stri += str(self.mrBayes_fn03_cov_var.get())  + '\t\n'
        stri += str(self.mrBayes_fn04_cov_var.get())  + '\t\n'        

        print stri

        filename = self.rootTable + self.organism  + '_' + self.gene_title + '_default.ini'
        self.write(filename, stri)

        if self.new_search:
            self.new_search = False
        else:
            self.filter_listSpecies(mode=self.rb_seq_filter_var.get(), cutoff=self.cutoffNumSeq_var.get())

    def call_external_program(self, program, param):
        call([program, param])
        
    def write(self, filename, stri):
                
        try:
            f = open(filename, 'w')
            f.write(stri)
            f.flush()
            
            # tkMessageBox.showinfo('message', 'saving: %s \n\nin %s'%(stri, self.filename))
            ok = True
        except:
            ok = False
            
        finally:
            if not ok:
                print 'Could not write in %s'%(filename)
            else:
                print 'Saved %s'%(filename)
                f.close()
            
        return ok
        
        
            
    def read_text_to_species_list(self, filename):
        ret, lines = self.read_text(filename)
        
        if not ret:
            return None
        
        return lines
        
        
        
    def read_text(self, filename):
        try:
            if not os.path.isfile(filename):
                print 'Could not find %s'%(filename)
                return False, None
            
            f = open(filename, 'r')
            
            lines = f.readlines()
            for i in range(len(lines)):
                lines[i] = lines[i].rstrip('\n')
            
            fileError = False
        except:
            fileError = True
            
         
        if fileError:
            print 'Could not find %s'%(filename)
            return False, None
                
        return True, lines
    
    def recalc(self, seqs):

        dic_species = {}
        
        for seq in seqs:
            desc = seq.id
            species = self.find_seq_species(desc)
            
            if species in dic_species.keys():
                dic_species[species].append(len(seq.seq))
            else:
                dic_species[species] = [len(seq.seq)]
                
        return dic_species

    def redefine_params(self, dic_seqs_mincut, dic_seqs_maxmer, which):
        iLoop = 0
        lista = self.speciesListbox.get(0,END)
    
        stri = ""
        for line in lista:
            iLoop += 1
            species = self.find_species(line)
            
            mat = self.dicParams[species]
            
            
            matmin = dic_seqs_mincut[species] if species in dic_seqs_mincut else None
            matmax = dic_seqs_maxmer[species] if species in dic_seqs_maxmer else None
            
            nSeq = mat[1]

            mat[5] = matmin[0] if matmin else 0
            mat[6] = matmax[0] if matmax else 0
            
            ''' at alignment look at original n (nSeq), in purge and consensus not '''
            if matmin and matmax:
                if which == 'align':
                    if nSeq == len(matmin) and nSeq == len(matmax):
                        continue
                else:
                    if mat[2] == len(matmin) and mat[3] == len(matmax):
                        continue
                
            mat[2] = len(matmin) if matmin else 0
            mat[3] = len(matmax) if matmax else 0
            
            stri += "%s changed from n=%i to n(mincut)=%i and n(maxmer)=%i\n"%(species, mat[1], mat[2], mat[3])
    
        self.filter_listSpecies(mode=self.rb_seq_filter_var.get(), cutoff=self.cutoffNumSeq_var.get())
        self.saveSpeciesListParams()
        
        if stri == '':
            stri = "There are no corrections."
            tkMessageBox.showinfo("All right", stri)
        else:
            tkMessageBox.showinfo("Warning", stri)
            
        print("")
        print(stri)
        print("")
 
    def find_species(self, desc):
        pos = desc.find("|")
        if pos < 1:
            return ""
        desc = desc[pos+1:]
        
        pos = desc.find("|")
        if pos < 1:
            return ""
        
        return desc[pos+1:].strip()
         
    def find_seq_species(self, desc):
        pos = desc.find("|")
        if pos < 1:
            return ""
        desc = desc[pos+1:]
        
        pos = desc.find("|")
        if pos < 1:
            return ""
        
        return desc[:pos].strip()
             
    def minimize_window(self, w=tk_root):
        w.iconify() # wm_state('iconic')
        
    def normalize_window(self, w=tk_root):
        w.deiconify() # wm_state('iconic')
        
class IODirector(object):
    def __init__(self,text_area):
        self.text_area = text_area

class StdoutDirector(IODirector):
    def write(self,stri):
        self.text_area.insert(tk.END,stri)
        self.text_area.see(tk.END)
        
    def flush(self):
        pass
    
tk_root = tk.Tk()
width, height = tk_root.winfo_screenwidth()*.98, tk_root.winfo_screenheight()*.9

x_center = (tk_root.winfo_screenwidth() // 2) - (width // 2)
y_center = (tk_root.winfo_screenheight() // 2) - (height // 2)

tk_root.geometry("%dx%d+%d+%d" % (width, height, x_center, y_center))

Desktop()
tk_root.mainloop()
