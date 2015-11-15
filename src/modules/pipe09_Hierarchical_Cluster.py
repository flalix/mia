#-*- coding: utf-8 -*-
'''
Created 08/08/2012
Updated 01/09/2014
Updated 11/10/2015 - get_params()

@author: Damian Kao
@url: http://blog.nextgenetics.net/?e=44

'''
import FileDialog
from scipy.cluster.hierarchy import linkage 
from scipy.cluster.hierarchy import dendrogram
import classes.BarGraphic as graphPac
import matplotlib.pyplot as plt
import classes.Drosophila as dro
import os
import numpy as np

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
        
        if desk.cluster_method == 'complete':
            desk.cluster_method_desc = desk.cluster_method + '-max'
        elif desk.cluster_method == 'single':
            desk.cluster_method_desc = desk.cluster_method + '-min'
        elif desk.cluster_method == 'weighted':
            desk.cluster_method_desc = desk.cluster_method 
        else:
            desk.cluster_method = 'centroid'
            desk.cluster_method_desc = desk.cluster_method 
            
        desk.cluster_method_desc = desk.cluster_method_desc[0].upper() + desk.cluster_method_desc[1:]

        if desk.organism == '':
            self.error_msg = "Define the organism or write 'any'"
            return

        if desk.gene == '' and desk.title == '':
            self.error_msg = 'Define at least one Gene or Title'
            return
        
        ''' -----------------------------------------'''
        if desk.each_all == 'each':
            matLoopList = [ [desk.withCorrection, desk.minmax]]
        else:
            matLoopList = [ [kCorr, kMinMax] for kCorr in [False, True] for kMinMax in ['mincut','maxmer']]

    
        for matLoop in matLoopList:
            if not self.looping(desk, matLoop):
                return

        self.failure = False
        self.error_msg = 'Task ended. All right.'
        return

    def looping(self, desk, opt):
        
        plt.close("all")
        plt.clf()
                    
        desk.withCorrection = opt[0]
        desk.minmax = opt[1]

        if desk.withCorrection:
            self.str_correction = '-bias corr.'
            self.filename_correction = '_bias_corr'
        else:
            self.str_correction = ''
            self.filename_correction = ''
    
        print("\n--->>>", desk.minmax, self.str_correction)
           
        desk.colorThreshold = desk.colorThreshold_var.get()
        
        if desk.mnat:
            desk.unit = 'mnat'
            desk.factor = 1000
            desk.roundVal = 2
        else:
            desk.unit = 'nat'        
            desk.factor = 1
            desk.roundVal = 4        

        if desk.vert_horiz == 'HMI':
            xLabel = 'JSD(HMI) (%s)'%(desk.unit)

            title = "Hierarchical Cluster Method=%s of JSD(HMI)- %s %s %s"\
             %(desk.cluster_method_desc, desk.organism, desk.seqType, desk.gene)
            
            if desk.frame > 0:
                title += '\nJSD(HMI) %s%s, desk.frame %i, #letter %i, min(L)=%i, min(#seqs)=%i' % \
                (desk.minmax, self.str_correction, desk.frame, desk.numOfLetters, desk.cutoffLength, desk.cutoffNumSeq)
            else:
                title += '\n%s%s, letter %i, min(L)=%i, min(#seqs)=%i' % \
                (desk.minmax, self.str_correction, desk.numOfLetters, desk.cutoffLength, desk.cutoffNumSeq)
        elif desk.vert_horiz == 'VMI':
            xLabel = 'JSD(VMI) (%s)'%(desk.unit)
                

            ''' multidimensional distance '''
            title = "Hierarchical Cluster Method=%s of JSD(VMI), %s %s %s"\
                %(desk.cluster_method_desc, desk.organism, desk.seqType, desk.gene)

            title += '\n%s%s, #letter %i, min(L)=%i, min(#seqs)=%i' % \
                (desk.minmax, self.str_correction, desk.numOfLetters, desk.cutoffLength, desk.cutoffNumSeq)

        else:
            xLabel = 'JSD(VSH) (nat)'

            ''' multidimensional distance '''
            title = "Hierarchical Cluster Method=%s of JSD(VSH), %s %s %s"\
                %(desk.cluster_method_desc, desk.organism, desk.seqType, desk.gene)

            title += '\n%s%s, #letter %i, min(L)=%i, min(#seqs)=%i' % \
                (desk.minmax, self.str_correction, desk.numOfLetters, desk.cutoffLength, desk.cutoffNumSeq)
        
        desk.set_cluster_filenames()
        filename = desk.cluster_input_filename

        
        ret, _, colHeaders, dataMatrix = self.open_distance_matrix_file(desk.rootTable + filename)
        if not ret:
            self.error_msg = 'Could not find %s'%(desk.rootTable + filename)
            return False
        
        pictureName = 'Cluster_' + filename.replace('.txt','')

        ''' desk.dr defined in pipe_desktop get_params() '''
        if desk.dr:
            rows = desk.dr.labels(colHeaders)
        else:
            rows = colHeaders
         
        #convert native python array into a numpy array
        # dataMatrix = log10(dataMatrix)
        # print dataMatrix

        dataMatrix = np.array(dataMatrix)
        maxDist = 0
        if desk.factor != 1:
            for i in range(len(dataMatrix)):
                for j in range(len(dataMatrix[i])):
                    dataMatrix[i][j] = dataMatrix[i][j] * desk.factor
                    if dataMatrix[i][j] > maxDist:
                        maxDist = dataMatrix[i][j]
        else:
            for i in range(len(dataMatrix)):
                for j in range(len(dataMatrix[i])):
                    if dataMatrix[i][j] > maxDist:
                        maxDist = dataMatrix[i][j]                    

        # single, weighted, average, co    mplete
        linkageMatrix = linkage(dataMatrix, method=desk.cluster_method, metric='euclidean')
        
        ''' finding maximum '''
        maxLinkDist = 0
        for i in range(len(linkageMatrix)):
            for j in range(len(linkageMatrix[i])):
                if linkageMatrix[i][j] > maxLinkDist:
                    maxLinkDist = linkageMatrix[i][j]  
        
        ''' hierarchical cluster distorce distances
        factor = maxDist/(2*maxLinkDist) '''
        
        for i in range(len(linkageMatrix)):
            linkageMatrix[i][2] = round(linkageMatrix[i][2]*.5, desk.roundVal)

                    
        fig = plt.figure(1, dpi=desk.dpi)
        ax = fig.add_subplot('111')

        plt.subplots_adjust(bottom=.1, left=.05, right=.84)
             
        yLabel = 'species'
        
        plt.rcParams['lines.linewidth'] = 2.5
        fontsize = 26
        
        plt.title(title, fontsize=fontsize) 
        ax.set_xlabel(xLabel, fontsize=fontsize)
        ax.set_ylabel(yLabel, fontsize=fontsize)

        
        # make colorbar labels bigger
        leaf_font_size = 28

                    
        ''' ddata = '''
        try:
            dendrogram(linkageMatrix,
                       color_threshold=desk.colorThreshold,
                       labels=rows, orientation='right')  # show_leaf_counts=True   , leaf_font_size=leaf_font_size
        except:
            print("Failed in printing dendrogram")
            pass

        plt.xticks(fontsize=leaf_font_size)
        plt.yticks(fontsize=leaf_font_size)
        
        '''
        # print ddata
        spList = ddata['ivl']
        # print len(spList), spList
        nickList = copy.deepcopy(spList)
        nickList.sort()
        
        dic = {}
        for i in range(len(spList)):
            sp = spList[i]
            for j in range(len(nickList)):
                if sp == nickList[j]:
                    dic[i] = j
                    #print i, spList[i], ' equal ',j, nickList[j]
                    break
        
        
        count = 0
        for i, d in zip(ddata['icoord'], ddata['dcoord']):
            count += 1
            # print i, d
        
            # specie01 x specie02 - mean error distance 
            num  = (i[0]-5)/10.
            sp1a = int(num)
            diff = num - sp1a
        
            if diff == 0:
                wei1a = 1
                sp1b  = sp1a
                wei1b = 0
            else:
                sp1b = sp1a+1
                wei1a = diff
                wei1b = 1. - wei1a
                
            #if num == 0:
            #    print '>>>> viri'
            num  = (i[2]-5)/10.
            sp2a = int(num)
            diff = num - sp2a
        
            if diff == 0:
                sp2b  = sp2a
                wei2a = 1
                wei2b = 0
            else:
                sp2b = sp2a+1
                wei2a = diff
                wei2b = 1. - wei2a
                
            #print sp1a, sp1b, sp2a, sp2b
            #print wei1a, wei1b, wei2a, wei2b
        
        
            ste = 0.
            if wei1a>0 and wei2a>0:
                ste += wei1a*wei2a*seMatrix[dic[sp1a]][dic[sp2a]]  
            if wei1a>0 and wei2b>0:
                ste += wei1a*wei2b*seMatrix[dic[sp1a]][dic[sp2b]]  
            if wei1b>0 and wei2a>0:
                ste += wei1b*wei2a*seMatrix[dic[sp1b]][dic[sp2a]]  
            if wei1b>0 and wei2b>0:
                # print sp1b, sp2b
                ste += wei1b*wei2b*seMatrix[dic[sp1b]][dic[sp2b]]
            
            ste = round(ste,4)
        
            dist = seMatrix[dic[sp1a]][dic[sp2a]]
            dist = round(dist,4)
            
            # print 'dist', dist, 'ste', ste
        
            x = 0.5 * sum(i[1:3])
            y = round(d[1],4)
            stry = str(y) + '\nd='+str(dist) + '\nse='+str(ste)
            plt.plot(x, y, 'ro')
            stry = ''
            if abs(y) > desk.colorThreshold:
                plt.annotate(stry, (x, y), xytext=(0, -8),
                             textcoords='offset points',
                             va='top', ha='center')
        '''
        
        self.myPlot = graphPac.Plot()
        self.myPlot.print_graph(desk, fig, pictureName=pictureName, frame=desk.tk_root, stay=True)
                    
        return True

    def open_distance_matrix_file(self, filename):
            #open the file assuming the data above is in a file called 'dataFile'
            if not os.path.exists(filename):
                return False, None, None, None
                
            # print 'opening %s'%(filename)
        
            try:
                inFile = open(filename,'r')
                #save the column/row headers (conditions/genes) into an array
                colHeaders = inFile.next().strip().split()[1:]
                rowHeaders = []
                dataMatrix = []
                
                for line in inFile:
                    data = line.strip().split('\t')
                    if data[0] == '':
                        break
                    rowHeaders.append(data[0])
                    dataMatrix.append([np.double(x) for x in data[1:]])
            except:
                return False, None, None, None
        
            return True, rowHeaders, colHeaders, dataMatrix
        
    