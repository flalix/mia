'''
Created on 22/01/2015

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica
'''
from Tkinter import END
from Bio import SeqIO, AlignIO
from subprocess import Popen, PIPE
import os, copy
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
        
class Pipe():
    def __init__(self, desk):
        self.desk = desk
        self.failure = True
        self.error_msg = ''

        try:
            desk.get_params()
            desk.showGraph = False
        except:
            self.error_msg = 'Could not get parameters.'
            return
        
     
        if desk.organism == '':
            self.error_msg = "Define the organism or write 'any'"
            return
        
        if desk.gene_title == '':
            self.error_msg = 'Define at least one Gene or Title'
            return
        
        
        if not desk.dicParams.keys():
            self.error_msg =  "Could not find list to perform the alignment."
            return

        ''' second loop - merging all sequences from all species ''' 
        allFilename = desk.rootFasta + ('%s_all_%s_%s_%iL_cutoff%i.fasta') %\
             (desk.organism, desk.seqType, desk.gene_title, desk.cutoffLength, desk.cutoffNumSeq)
        allFilenameAligned = desk.rootFasta + ('%s_all_%s_%s_%iL_cutoff%i_aligned.fasta') %\
             (desk.organism, desk.seqType, desk.gene_title, desk.cutoffLength, desk.cutoffNumSeq)


        ''' first loop - aligning each species ''' 
        iLoop = 0
        lista = desk.speciesListbox.get(0,END)
    
        for line in lista:
            iLoop += 1
            species = desk.find_species(line)

            mat = desk.dicParams[species]
            if mat[0] != 'x':
                desk.showmsg_obs('%i/%i - %s not chosen - not checked.'%(iLoop, len(lista), species))
                continue
            
            inFile = desk.rootFasta + ('%s_%s_%s_%s_%iL') %\
                 (desk.organism, species, desk.seqType, desk.gene_title, desk.cutoffLength)
                             
            outFile = inFile + '_aligned.fasta'
            inFile += '.fasta'
            
            if os.path.exists(outFile) and not desk.realing_each_species:
                desk.showmsg_obs("Muscle Alignment for %s already exists."%(outFile))
            else:
                desk.showmsg_obs("Muscle Alignment for %s %s"%(desk.organism, species))
    
                seqs_align = self.align_muscle_dna(desk, inFile)
                
                if not seqs_align:
                    self.error_msg = 'Error while reading/aligning %s.'%(inFile)
                    return
                
                if not self.save_fasta(outFile, seqs_align):
                    self.error_msg = 'Error while writing %s.'%(inFile)
                    return
                
                        
        desk.showmsg_obs("Second loop, merging all sequences.")

        if not os.path.exists(allFilenameAligned) or desk.realing_all:

            ''' second loop - merging all sequences from all species ''' 
            desk.showmsg_obs("Merging all aligned sequences.")
            all_seqs = []

            iLoop = 0
            
            ''' reading each species '''
            for line in lista:
                iLoop += 1
                species = desk.find_species(line)
                
                mat = desk.dicParams[species]
                if mat[0] != 'x':
                    desk.showmsg_obs('%i/%i - %s not chosen - not checked.'%(iLoop, len(lista), species))
                    continue                
    
                inFile = desk.rootFasta + ('%s_%s_%s_%s_%iL_aligned.fasta') %\
                     (desk.organism, species, desk.seqType, desk.gene_title, desk.cutoffLength)
                
                # desk.showmsg_obs("%s len=%i"%(inFile, numOfSeqs))
                
                seqs = self.read_fasta(inFile)
                
                if not seqs:
                    self.error_msg = 'Error while reading %s.'%(inFile)
                    return
                
                ''' merging in all_seqs '''
                for seq in seqs:
                    all_seqs.append(seq)
    
            desk.showmsg_obs("Saving all merged sequences.")
            if not self.save_fasta(allFilename, all_seqs):
                self.error_msg = 'Error while writing %s.'%(inFile)
                return
            
            desk.showmsg_obs("Aligning all sequences with muscle.")

            seqs_align = self.align_muscle_dna(desk, allFilename)
            
            if not seqs_align:
                self.error_msg = 'Error while reading/aligning all sequences: %s.'%(inFile)
                return
            
            desk.showmsg_obs("Saving: %s"%allFilenameAligned)
            if not self.save_fasta(allFilenameAligned, seqs_align):
                self.error_msg = 'Error while writing all sequences %s.'%(allFilenameAligned)
                return
        
         
        '''
            creating mincut - maxmer
        '''
        self.mincutFile = desk.rootFasta + ('%s_mincut_%s_%s_%iL_cutoff%i_aligned.fasta') %\
                 (desk.organism, desk.seqType, desk.gene_title, desk.cutoffLength, desk.cutoffNumSeq)
        self.maxmerFile = desk.rootFasta + ('%s_maxmer_%s_%s_%iL_cutoff%i_aligned.fasta') %\
                 (desk.organism, desk.seqType, desk.gene_title, desk.cutoffLength, desk.cutoffNumSeq)
                 
        '''    maxmer with minimum filter '''
        desk.showmsg_obs('Creating maxmer file with Vertical minimum cutoff.')
        seqs = self.read_fasta(allFilenameAligned)
        
        if not seqs:
            self.error_msg = 'Error while reading %s.'%(allFilenameAligned)
            return
                             
        seqs2 = copy.deepcopy(seqs)
        seqs = self.cut_vertical_gaps(seqs, "maxmer")
        seqs = self.cut_horizotal_gaps(seqs)
        
        if not self.save_fasta(self.maxmerFile, seqs):
            self.error_msg = 'Error while writing mincut file %s.'%(self.maxmerFile)
            return

        dic_seqs_maxmer = desk.recalc(seqs)

        '''    mincut with maximum filter '''
        desk.showmsg_obs('Creating mincut file with Vertical maximum cutoff.')
        
                      
        seqs = self.cut_vertical_gaps(seqs2, "mincut")
        seqs = self.cut_horizotal_gaps(seqs)

        if not self.save_fasta(self.mincutFile, seqs):
            self.error_msg = 'Error while writing maxmer file %s.'%(self.mincutFile)
            return
        
        dic_seqs_mincut = desk.recalc(seqs)

        print("Aligned fasta files created. Verify sequences with an editor, and go to Purging")
        print("all sequeces file: %s"%allFilenameAligned)
        print("mincut file: %s"%self.mincutFile)
        print("all maxmer file: %s"%self.maxmerFile)
        desk.showmsg_obs("")
        
        desk.redefine_params(dic_seqs_mincut, dic_seqs_maxmer, 'align')
        self.failure = False


    def align_muscle_dna(self, desk, inFile):
        '''
        muscle -in input file (fasta) [-out output file (default fasta)]
              [-diags] [-log log file] [-maxiters n] [-maxhours n] [-maxmb m]
              [-html] [-msf] [-clw] [-clwstrict] [-log[a] logfile] [-quiet]
              [-stable] [-group] [-version]
              
       -stable
           Output sequences in input order (default is -group)

       -group
           Group sequences by similarity (this is the default)
        '''

        if desk.alignment_var.get() == 'Muscle':
            cmd = ['muscle', "-quiet", "-maxiters", "1", "-diags"]
        else:
            cmd = ['clustal']
            
        if desk.isWindows:

            muscle_exe = r"C:\seaview\muscle.exe" 
            if not os.path.isfile(muscle_exe):
                print("%s executable missing"%(muscle_exe))
                return None
                        
            muscle_cline = MuscleCommandline(muscle_exe, input=inFile.replace("/","\\"))
            try:
                # stdout, stderr = muscle_cline()
                stdout, _ = muscle_cline()
            except:
                print("\n---------------------------")
                print("Error while aligning with Muscle. Do it manually.")
                print("---------------------------\n")
                return None
            
            seqs_align = AlignIO.read(StringIO(stdout), "fasta")
                
        
            '''
            try:
                # muscle_cline = MuscleCommandline(muscle_exe, input=inFile, out=outFile)
                muscle_cline = MuscleCommandline(muscle_exe, input=inFile.replace("/","\\"))
                stdout, stderr = muscle_cline()
                seqs_align = AlignIO.read(StringIO(stdout), "fasta")
            except:
                return []
            '''
        else:
            muscle = Popen(cmd, stdin=PIPE, stdout=PIPE, universal_newlines=True)

            seqs = self.read_fasta(inFile)
            
            if not seqs:
                print 'Error while reading fasta'
                return None
            
            SeqIO.write(seqs, muscle.stdin, "fasta")  # Send sequences to Muscle in FASTA format.
            muscle.stdin.close()
    
            seqs_align = AlignIO.read(muscle.stdout, 'fasta')  # Capture output from muscle and get it into FASTA format in an object.
    
            muscle.stdout.close()

        return seqs_align


    def cut_vertical_gaps(self, seqs_align, which):
        self.desk.showmsg_obs("Alignment done, cutting vertical gaps")
        
        if which == "mincut":
            cutoff = self.desk.minVertCutoff/100.
        else:
            cutoff = self.desk.maxVertCutoff/100.
        
        if cutoff == 0:
            return seqs_align
        
        seqs = []
        for seq in seqs_align:
            seqs.append(seq)
            
        seqLen = len(seqs[0])
        numSeqs = len(seqs)
        
        maxGap = numSeqs * cutoff
        gaps = []
                    
        for col in range(seqLen):
            count = sum([1 for row in range(numSeqs) if seqs[row].seq[col] == '-'])

            if count > maxGap:
                gaps.append(col)

        if not gaps:
            return seqs_align
        
        this = -1
        replaces = []
        
        for k in gaps:
            this += 1
            
            if k == this:
                continue
        
            replaces.append([this,k])
            this = k
        
        this = gaps[len(gaps)-1] + 1
        if this < seqLen:
            replaces.append([this,seqLen])
            
        
        # print 'Selected bp regions', replaces
    
            
        for row in range(numSeqs):
            seq = seqs[row]
            seq2 = ''
            for limInf, limSup in replaces:
                seq2 += seq.seq[limInf: limSup]
                
            seqs_align[row].seq = seq2
        
        return seqs_align
        
    def cut_horizotal_gaps(self, seqs_align):
        self.desk.showmsg_obs("Cutting horizontal gaps")
        cutoff = self.desk.horizCutoff/100.
        if cutoff == 0:
            return seqs_align

        seqLen = len(seqs_align[0])
        numSeqs = len(seqs_align)
        
        maxGap = seqLen * cutoff
        
        badList = [ i for i in range(numSeqs) if  seqs_align[i].seq.count('-') > maxGap]

        '''
        badList = []
        for i in range(numSeqs):
            count = seqs_align[i].seq.count('-')
                    
            if count > maxGap:
                badList.append(i)
        '''
        
        offset = 0
        for k in badList:
            seqs_align.pop(k-offset)
            offset += 1
        
        return seqs_align
        
       
                
    def read_fasta(self, inFile):
        seqs = []
            
        try:
            for seq in SeqIO.parse(inFile, format="fasta"):
                seqs.append(seq)
        except:
            print 'coultd not read %s'%(inFile)


        return seqs
    
    def save_fasta(self, filename, seqs):
        
        try:
            handle = open(filename, "w")
            SeqIO.write(seqs, handle, 'fasta')
            ret = True
            print 'Saved', filename
            
        except:
            print "File '%s' not saved. Writing error. See if seq_record is well formated."%(filename)
            ret = False

        finally:   
            handle.close()
            
        return ret

    
    '''
    def align_muscle2(self, inFile, outFile):
        from Bio.Align.Applications import MuscleCommandline
        muscle_cline = MuscleCommandline(input=inFile, out=outFile, "-stable")
        stdout, stderr = muscle_cline()
        align = AlignIO.read(StringIO(stdout), "fasta")
        print(align)
        
        

        
    def align_smith_waterman(self, inFile, outFile):
        # Smith-Waterman algorithm local alignment, and Needleman-Wunsch global alignment
        from Bio.Emboss.Applications import WaterCommandline
        water_cmd = WaterCommandline(gapopen=10, gapextend=0.5,\
                                     stdout=True, auto=True,\
                                     asequence="a.fasta", bsequence="b.fasta")
        print("About to run: %s" % water_cmd)
        std_output, err_output = water_cmd()
        
        return std_output, err_output
    

        from Bio.Align.Applications import MuscleCommandline
        from StringIO import StringIO

        muscle_exe = r"C:\seaview\muscle.exe" 
        if not os.path.isfile(muscle_exe):
            print("%s executable missing"%(muscle_exe))
            return []

    
        try:
            # muscle_cline = MuscleCommandline(muscle_exe, input=inFile, out=outFile)
            muscle_cline = MuscleCommandline(muscle_exe, input=inFile.replace("/","\\"))
            stdout, stderr = muscle_cline()
            seqs_align = AlignIO.read(StringIO(stdout), "fasta")
        except:
            return []
                
                    
        if Clustalw_Muscle == 'C':        
            from Bio.Align.Applications import ClustalwCommandline
            
            print 'filename', mSeq.filename
        
            clustalw_exe = r"C:\python\seaview4\clustalw2.exe"
            print 'generates automatically aln file',mSeq.fileNameAln,'from', mSeq.completeName
            clustalw_cline = ClustalwCommandline(clustalw_exe, infile=mSeq.completeName)
            assert os.path.isfile(clustalw_exe), "Clustal2.exe executable missing"
            stdout, stderr = clustalw_cline()
            
            try:
                align = AlignIO.read(mSeq.fileNameAln, "clustal")
                print 'clustal alignment'
                print align    
            
                tree = Phylo.read(mSeq.fileNameDnd, format="newick")
                Phylo.draw_ascii(tree)
            except:
                print "couldn't read", mSeq.fileNameAln
                
            
        else:
                
        from Bio.Align.Applications import MuscleCommandline
        
        muscle_exe = "C:\python\seaview4\muscle.exe" 
        in_file = mSeq.completeName
        out_file = mSeq.fileNameMscAln
        print 'muscle: generates automatically aln file',mSeq.fileNameAln,'from', mSeq.completeName
        muscle_cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file) 
    
        assert os.path.isfile(muscle_exe), "Muscle.exe executable missing"
        stdout, stderr = muscle_cline()
        
        # os.rename(mSeq.filename, mSeq.filename + '_ori')
        try:
            align = AlignIO.read(handle=out_file, format="fasta")
            print 'muscle alignment'
            print align
    
        except:
            print "couldn't read", out_file
            exit() 
        
'''