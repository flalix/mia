#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on Set 09, 2011
Updated on May 13, 2014
Updated on Oct 02, 2015

@author: Flavio Lichtenstein
'''
# no module named FileDialog  # seemed to be required also by matplotlib
# and after I finally added import FileDialog the program actually works!
# import FileDialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d import Axes3D
# import importlib
# importlib.import_module('mpl_toolkits.mplot3d').__path__

import numpy as np
import Tkinter as tk
from scipy.cluster.hierarchy import linkage, dendrogram

class Plot:
    def __init__(self):
        pass
    
    def print_graph(self, desk, fig, pictureName, frame=None, stay=False):
        pictureName = desk.rootImage + pictureName + '.' + desk.imgType
        # fig.tight_layout()
        
        if desk.isWindows or not frame:
            self.plot_plot(desk, pictureName, stay=stay)
        else:
            self.plot_tkinter(desk, fig, frame, pictureName)
            


    def plot_plot(self, desk, pictureName, stay=False):
        mng = plt.get_current_fig_manager()
        '''
        if desk.isWindows:
            mng.window.state('zoomed') #works fine on Windows!
        '''

        if plt.get_backend() == 'TkAgg': 
            try:
                mng.window.state('zoomed')  # #works fine on Windows!
            except:
                try:
                    mng.frame.Maximize(True)
                except:
                    try:
                        mng.window.showMaximized() 
                    except:
                        # print "Could not maximize"
                        mng.window.wm_geometry("1400x900+50+50")
       
        if stay and not desk.saveGraph:                       
            plt.show()
        else:
            plt.show(block=False)
            
        if desk.saveGraph:
            plt.savefig(pictureName, format=desk.imgType, dpi=desk.dpi)

        if not desk.showGraph:
            plt.close()

        
    def savePlot(self, desk, pictureName):
        print("Saving image(%s): %s"%(desk.imgType, pictureName))
        plt.savefig(pictureName, format=desk.imgType, dpi=desk.dpi)  # , bbox_inches='tight'


    def plot_tkinter(self, desk, fig, frame, pictureName):
        '''
        if desk.saveGraph and not desk.showGraph:
            plt.draw()      # force a draw
            plt.savefig(pictureName, format=desk.imgType, dpi=desk.dpi)

        if not desk.showGraph:
            return
        '''
        
        master = tk.Toplevel(frame)
              
        '''
        if not desk.showGraph and desk.saveGraph:
            desk.minimize_window(master)      
        '''
        '''
            root.update()
            
            print (root.winfo_width())
            print (root.winfo_height())
            print (root.winfo_geometry())        
        '''
        # frame.update()
        master.geometry("%dx%d%+d%+d" % (frame.winfo_width(), frame.winfo_height(), 0, 0))

        '''
        master.geometry("%dx%d%+d%+d" % (1300, 700, 0, 0))

        master.columnconfigure(0, weight=1)
        master.rowconfigure(width=1, height=1)      
        master.wm_state('zoomed')
        aster.wm_title(self.title)
        '''

        canvas = FigureCanvasTkAgg(figure=fig,master=master)
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        ''' always show, otherwise saves badly, does not expand the figure ... 
        if desk.showGraph:
            canvas.show()
        '''
        canvas.show()

        #from time import sleep
        #sleep(0.1)
        
        if desk.saveGraph:
            self.savePlot(desk, pictureName)

        if not desk.showGraph:
            master.withdraw()
        
        '''
        del master
        '''
    
class BarGraphic:
    def __init__(self, dpi=300):
        self.dpi = dpi
        self.fontsize = 16- 2* int(round( (dpi-100)/ 100.,0))
        

    def gHist(self, lines, columns, numOfFig, havexSeq, xSeq, ySeq, yLine, parSDV, maxEntropy, par_xmin, par_xmax, par_yMin, par_ymax, par_xlabel, par_ylabel, par_title, par_width, par_color, wantTicks=True, seqSdv=None):
        plt.subplot(lines, columns, numOfFig)
        
        plt.subplots_adjust(left=0.10, right=0.95, bottom=0.05, top=0.95)
        
        wid = par_width / 35
        
        if (len(yLine) == 0):
            if (parSDV == 0):
                plt.bar(xSeq, ySeq, width=wid)
            else:
                plt.bar(xSeq, ySeq, yerr=parSDV, width=wid, error_kw=dict(elinewidth=wid / 3, color=par_color, ecolor='red'))
        else:
            plt.bar(xSeq, ySeq, yerr=parSDV, width=wid, error_kw=dict(elinewidth=wid / 3, color=par_color, ecolor='red'))  
            plt.plot(xSeq, yLine, color='red')
            
            seqSup = []
            seqInf = []
            seqMean = []
            
            for _ in range(len(yLine)):
                seqSup.append(maxEntropy)
                seqMean.append(maxEntropy - parSDV)
                seqInf.append(maxEntropy - 2 * parSDV)
                
            plt.plot(xSeq, seqSup, color='black')
            plt.plot(xSeq, seqMean, color='red')
            plt.plot(xSeq, seqInf, color='black')

        if wantTicks and havexSeq:
            plt.xticks(xSeq)

        self.yMax = par_ymax
        self.yMin = par_yMin
    
        self.xMax = par_xmax
        self.xMin = par_xmin   
        
        plt.xlim(par_xmin, par_xmax)
        plt.ylim(par_yMin, par_ymax)
        
        plt.title(par_title, fontsize=self.fontsize)
        plt.xlabel(par_xlabel, fontsize=self.fontsize)
        plt.ylabel(par_ylabel, fontsize=self.fontsize)


    # typeGraph par_width, yLine,parSDV,  maxEntropy,par_xmax
    def gBar(self, lines, columns, numOfFig, havexSeq, xSeq, ySeq, par_xmin, par_yMin, par_ymax, par_xlabel, par_ylabel, par_title, wantTicks=True, seqSdv=None):
        plt.subplot(lines, columns, numOfFig)
        
        plt.subplots_adjust(left=0.10, right=0.95, bottom=0.05, top=0.95)

        par_xmax = np.max(xSeq)
        plt.plot(xSeq, ySeq, 'r--', color='blue')
        if seqSdv:
            plt.errorbar(x=xSeq, y=ySeq, yerr=seqSdv, ecolor='red')

        if wantTicks and havexSeq:
            plt.xticks(xSeq)

        plt.xlim(par_xmin, par_xmax)
        plt.ylim(par_yMin, par_ymax)

        plt.title(par_title, fontsize=self.fontsize)
        plt.xlabel(par_xlabel, fontsize=self.fontsize)
        plt.ylabel(par_ylabel, fontsize=self.fontsize)


    def sameBar(self, lines, columns, numOfFig, xSeq, ySeq, linestyleCode, color):
        plt.subplot(lines, columns, numOfFig)

        yMax = np.max(ySeq)
        yMin = np.min(ySeq)
        
        if yMax > self.yMax:
            self.yMax = yMax
            
            if self.yMax >= 0:
                self.yMax *= 1.1
            else:
                self.yMax *= .9
                
        if yMin < self.yMin:
            self.yMin = yMin
            
            if self.yMin >= 0:
                self.yMin *= .9
            else:
                self.yMin *= 1.1
                            
        plt.ylim(self.yMin, self.yMax)        
        plt.plot(xSeq, ySeq, linestyleCode, color=color) 
            
    def gNestHist(self, lines, columns, numOfImage, xSeq, ySeq, yLine, parError, par_xmin, par_xmax, par_yMin, par_ymax, par_xlabel, par_ylabel, par_title, par_width, par_color, meanF=None, stdF=None, medianF=None):
        ax = plt.subplot(lines, columns, numOfImage)
        
        if (len(yLine) == 0):
            ax.bar(xSeq, ySeq, yerr=parError, width=par_width, error_kw=dict(elinewidth=par_width/3, color=par_color, ecolor='red'))
        else:
            ax.bar(xSeq, ySeq, yerr=parError, width=par_width, error_kw=dict(elinewidth=par_width/3, color=par_color, ecolor='red'))  

            plt.plot(xSeq, yLine, color='red')
            seqSup = []
            seqInf = []
            for i in range(len(yLine)):
                seqSup.append(yLine[i] + parError)
                seqInf.append(yLine[i] - parError)
            plt.plot(xSeq, seqSup, color='black')
            plt.plot(xSeq, seqInf, color='black')
            
        # vertical line from (70,100) to (70, 250)
        if meanF:
            plt.plot([meanF, meanF], [par_yMin, par_ymax], 'k-', lw=2, color='black')
            plt.plot([meanF+stdF, meanF+stdF], [par_yMin, par_ymax], '--', lw=2, color='red')
            plt.plot([meanF-stdF, meanF-stdF], [par_yMin, par_ymax], '--', lw=2, color='red')
            plt.plot([medianF, medianF], [par_yMin, par_ymax], 'k-', lw=2, color='yellow')
            
            ax.annotate(r'$1\sigma$', xy=((meanF+stdF)*1.05, (par_ymax-par_yMin)*.75), color='red')

        if ((par_yMin < 0.5) and (par_yMin > 0)):
            par_yMin = 0

        plt.xlim(par_xmin, par_xmax)
        plt.ylim(par_yMin, par_ymax)
        plt.title(par_title, fontsize=self.fontsize)
        
        plt.xlabel(par_xlabel, fontsize=self.fontsize)
        plt.ylabel(par_ylabel, fontsize=self.fontsize)
        
    def gNextBar(self, lines, columns, numOfImage, xSeq, ySeq, par_xmin, par_xmax, par_yMin, par_ymax, par_xlabel, par_ylabel, par_title, par_width, par_color, meanF=None, stdF=None, medianF=None):
        plt.subplot(lines, columns, numOfImage)
        
        plt.plot(xSeq, ySeq, 'r--', color='blue')

        if ((par_yMin < 0.5) and (par_yMin > 0)):
            par_yMin = 0

        plt.xlim(par_xmin, par_xmax)
        plt.ylim(par_yMin, par_ymax)
        plt.title(par_title, fontsize=self.fontsize)
        
        plt.xlabel(par_xlabel, fontsize=self.fontsize)
        plt.ylabel(par_ylabel, fontsize=self.fontsize)
        
        
    def printBar(self):
        plt.show()



class MultiLineWithTitle:
    def __init__(self, numLines, numCols, left, top, legColumns, legendTitle='', title='',dpi=120):

        self.myPlot = Plot()
        
        self.fig = plt.figure(1, dpi=dpi)

        self.numLines = numLines
        self.numCols = numCols
        self.left = left
        self.top = top
        self.legColumns = legColumns 
        
        self.legendTitle = legendTitle
        ''' 100 dpi = fs 10, 200 = fs 9, 300 = fs 8 '''
        self.fontsize = 10- 2* int(round( (dpi-100)/ 100.,0))

        self.fig.text(.5, .95, title, ha='center', fontsize=self.fontsize, color="blue")

        self.colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k',
                       'aqua', 'teal', 'lightslategray',
                       'maroon', 'chartreuse', 'fuchsia',
                       'navy', 'purple', 'silver', 'violet',
                       'b', 'g', 'r', 'c', 'm', 'y', 'k']


    def multilineGeneralWithTitle(self, iLoop, numFigure, xSeq, meanCurve3val, title, label, printLegend=True, xlabel='', ylabel=''):
        which = str(self.numLines) + str(self.numCols) + str(numFigure)
        ax = self.fig.add_subplot(which)


        '''
        wspace = 0.2   # the amount of width reserved for blank space between subplots
        hspace = 0.2   # the amount of height reserved for white space between subplots
        '''
        plt.subplots_adjust(left=0.1, right=0.90, bottom=0.15, top=0.85, wspace = 0.4)
        
        ax.set_title(title, fontsize=self.fontsize)
        ax.set_xlabel(xlabel, alpha=0.5, fontsize=self.fontsize-1)
        ax.set_ylabel(ylabel, alpha=0.5, fontsize=self.fontsize-1)

        #ax.tick_params(axis='both', which='major', labelsize=self.fontsize-2)
        #ax.tick_params(axis='both', which='minor', labelsize=self.fontsize-2)
        plt.tick_params(labelsize=self.fontsize-2)

        try:
            color1 = self.colors[iLoop]
        except:
            color1 = 'blue'
        
        value = []
        sdvCurve = []
        for i in range(len(meanCurve3val)):
            x = meanCurve3val[i]
            value.append(x[0])
            sdvCurve.append(x[2])

        if (color1):
            ax.plot(xSeq, value, color=color1, label=label)
            # http://www.thetechrepo.com/main-articles/469-how-to-change-line-properties-in-matplotlib-python
            ax.plot(xSeq, sdvCurve, color=color1, linestyle='-')

        else:
            ax.plot(xSeq, value, label=label)
            ax.plot(xSeq, sdvCurve, linestyle='dotted')

        if printLegend:
            ax.legend(loc="center left", bbox_to_anchor=[self.left, self.top],
                       ncol=self.legColumns, shadow=True, title=self.legendTitle)
            ax.get_legend().get_title().set_color("blue")
        

        
    def multilineGeneral_TR_WithTitle(self, iLoop, numFigure, xSeq, entropyCurve3Val, title, label, printLegend=True, xlabel='', ylabel=''):
        which = str(self.numLines) + str(self.numCols) + str(numFigure)
        ax = self.fig.add_subplot(which)

        ax.set_title(title)
        ax.set_xlabel(xlabel, alpha=0.5, fontsize=self.fontsize)
        ax.set_ylabel(ylabel, alpha=0.5, fontsize=self.fontsize)
                  
        try:
            color1 = self.colors[iLoop]
        except:
            color1 = None
        
        meanCurveMax = []
        meanCurve = []

        # Tsallis and Renyi have errors
        for i in range(len(entropyCurve3Val)):
            x = entropyCurve3Val[i]
            print 'x', x
            value = x[0]
            sdv   = x[2]
            
            print 'value', value, 'sdv', sdv
            
            meanCurveMax.append(sdv)
            meanCurve.append(value)
                
        if (color1):
            ax.plot(xSeq, meanCurveMax, color=color1, linestyle='dotted')
            ax.plot(xSeq, meanCurve, color=color1, label=label)

        else:
            ax.plot(xSeq, meanCurveMax, linestyle='dotted')
            ax.plot(xSeq, meanCurve, label=label)

        if printLegend:
            ax.legend(loc="center left", bbox_to_anchor=[self.left, self.top],
                       ncol=self.legColumns, shadow=True, title=self.legendTitle)
            ax.get_legend().get_title().set_color("blue")
            
                    
    def multilineGeneralWithTitle_and_HorizError(self, iLoop, numFigure, xSeq, meanCurve, error, title, label, printLegend=True, xlabel='', ylabel=''):
        
        which = str(self.numLines) + str(self.numCols) + str(numFigure)
        ax = self.fig.add_subplot(which)

        ax.set_title(title)
        ax.set_xlabel(xlabel, alpha=0.5, fontsize=self.fontsize)
        ax.set_ylabel(ylabel, alpha=0.5, fontsize=self.fontsize)
        
        try:
            color1 = self.colors[iLoop]
        except:
            color1 = None
        
        if (color1):
            ax.plot(xSeq, meanCurve, color=color1, label=label)
        else:
            ax.plot(xSeq, meanCurve, label=label)

        if printLegend:
            ax.legend(loc="center left", bbox_to_anchor=[self.left, self.top],
                       ncol=self.legColumns, shadow=True, title=self.legendTitle)
            ax.get_legend().get_title().set_color("blue")
        
            
        for i in range(len(xSeq)):
            x = xSeq[i]
            y = meanCurve[i]
            size = error[i]
                
            square = plt.Rectangle((x - size / 2, y - size / 2), size, size, facecolor='gray', color=color1)
            ax.add_patch(square)
            
        horZ = []
        print 'xSeq', np.max(xSeq), '  error', np.max(error)
        num = round((np.max(xSeq) + np.max(error)) / 10., 4)
        val = -num
        for i in range(11):
            horZ.append(val)
            val += num
        plt.xticks(horZ)
        ax.grid(True)

                
    def multilineGeneralWithTitle_and_VertError(self, iLoop, numFigure, xSeq, incX, meanCurve, error, title, label, printLegend=True, xlabel='', ylabel='', goOrigem=False):
        
        which = str(self.numLines) + str(self.numCols) + str(numFigure)
        ax = self.fig.add_subplot(which)

        ax.set_title(title)
        ax.set_xlabel(xlabel, alpha=0.5, fontsize=self.fontsize)
        ax.set_ylabel(ylabel, alpha=0.5, fontsize=self.fontsize)
                  
        try:
            color1 = self.colors[iLoop]
        except:
            color1 = None
        
        if (color1):
            ax.plot(xSeq, meanCurve, color=color1, label=label)
        else:
            ax.plot(xSeq, meanCurve, label=label)

        if printLegend:
            ax.legend(loc="center left", bbox_to_anchor=[self.left, self.top],
                       ncol=self.legColumns, shadow=True, title=self.legendTitle)
            ax.get_legend().get_title().set_color("blue")
        
        
        horZ = []

        numPos = np.max(xSeq)
        numNeg = np.min(xSeq)
        
        if goOrigem:
            if numNeg < .5 and numNeg >= 0:
                numNeg = 0
            if numPos > -.5 and numPos < .5:
                numPos = 0    
        
        tot = numPos - numNeg
        widthX = tot / 1000. # 0.2 / (iLoop+1)
        maxi = np.max(meanCurve)
        withY = maxi / 1000.
                
        for i in range(len(xSeq)):
            x = xSeq[i]
            y = meanCurve[i]
            size = error[i]
            
            square = plt.Rectangle((x - widthX, y - size / 2.), 2.*widthX, size, facecolor='gray', color=color1)
            ax.add_patch(square)

            square = plt.Rectangle((x - (5 * widthX), y - size / 2.), (10.*widthX), withY, facecolor='gray', color=color1)
            ax.add_patch(square)

            square = plt.Rectangle((x - (5 * widthX), y + size / 2.), (10.*widthX), withY, facecolor='gray', color=color1)
            ax.add_patch(square)

        if incX:
            inc = incX
        else:
            inc = abs(round(tot / 20., 2))
        
        val = numNeg 
        for i in range(int(round(tot / inc, 0) + 1)):
            horZ.append(val)
            val += inc
            
        plt.xticks(horZ)
        ax.grid(True)

          
    
class MultiLine:
    def __init__(self, numLines, numCols, left, top, legColumns, legendTitle):
        self.numLines = numLines
        self.numCols = numCols
        self.left = left
        self.top = top
        self.legColumns = legColumns 
        
        self.legendTitle = legendTitle
        
        self.colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k',
                       'aqua', 'teal', 'lightslategray',
                       'maroon', 'chartreuse', 'fuchsia',
                       'navy', 'purple', 'silver', 'violet', 'white']
        '''
        b : blue
        g : green
        r : red
        c : cyan
        m : magenta
        y : yellow
        k : black
        w : white
        
        http://www.tedmontgomery.com/tutorial/colors.html

        '''
        
    ''' iLoop goes from 0 to numLines -1 '''
    def multiline(self, iLoop, numFigure, meanCurve, vc, prot, title):
        
        self.ax = plt.subplot(1, self.numCols, numFigure, title=title)
         
        try:
            color1 = self.colors[iLoop]
        except:
            color1 = None

        ''' http://matplotlib.sourceforge.net/api/pyplot_api.html '''
        label = prot.name.replace('DROSOPHILA ', '') + '(ind=' + str(prot.numIndiv) + ' vc=' + str(round(vc * 100, 3)) + '%)'
        
        if (color1):
            meanCurveAux = []
            if (vc > 0):
                for i in range(len(meanCurve)):
                    meanCurveAux.append(meanCurve[i] * (1 + vc))
                self.ax.plot(prot.alfa, meanCurveAux, color=color1, linestyle='dotted')

            self.ax.plot(prot.alfa, meanCurve, color=color1, label=label)

            if (vc > 0):
                for i in range(len(meanCurve)):
                    meanCurveAux[i] = meanCurve[i] * (1 - vc)
                self.ax.plot(prot.alfa, meanCurveAux, color=color1, linestyle='dotted')
        else:
            meanCurveAux = []
            for i in range(len(meanCurve)):
                meanCurveAux.append(meanCurve[i] * (1 + vc))
            self.ax.plot(prot.alfa, meanCurve, label=prot.name.replace('DROSOPHILA ', ''))
            self.ax.plot(prot.alfa, meanCurveAux, linestyle='dotted')
            for i in range(len(meanCurve)):
                meanCurveAux[i] = meanCurve[i] * (1 - vc)
            self.ax.plot(prot.alfa, meanCurveAux, linestyle='dotted')

        if (numFigure == 1):    
            self.ax.legend(loc="center left", bbox_to_anchor=[self.left, self.top],
                       ncol=self.legColumns, shadow=True, title=self.legendTitle)
            self.ax.get_legend().get_title().set_color("blue")
            
        self.ax.set_backgroundcolor(facecolor='b', alpha=0.5)
        plt.xticks(prot.alfa)
       
        # print 'iLoop', iLoop, 'numFigure', numFigure
        if (iLoop == (self.numLines - 1) and (numFigure == self.numCols)):   
            plt.draw()
            plt.show()       

                   
    def multilineGeneral(self, iLoop, numFigure, xSeq, meanCurve, vc, title, label, printLegend=True, xlable='', ylable=''):
        self.ax = plt.subplot(self.numLines, self.numCols, numFigure, title=title)
        plt.ylabel(ylable)
        plt.xlabel(xlable)

        try:
            color1 = self.colors[iLoop]
        except:
            color1 = None
        
        if (numFigure == 1):
            print 'color', color1, 'iLoop', iLoop, 'figure', numFigure, 'title', title

        ''' http://matplotlib.sourceforge.net/api/pyplot_api.html '''
        if (color1):
            meanCurveAux = []
            if (vc > 0):
                for i in range(len(meanCurve)):
                    meanCurveAux.append(meanCurve[i] * (1 + vc))
                self.ax.plot(xSeq, meanCurveAux, color=color1, linestyle='dotted')

            # print 'mandei bala', xSeq, meanCurve
            self.ax.plot(xSeq, meanCurve, color=color1, label=label)

            if (vc > 0):
                for i in range(len(meanCurve)):
                    meanCurveAux[i] = meanCurve[i] * (1 - vc)
                self.ax.plot(xSeq, meanCurveAux, color=color1, linestyle='dotted')
        else:
            meanCurveAux = []
            if (vc > 0):
                for i in range(len(meanCurve)):
                    meanCurveAux.append(meanCurve[i] * (1 + vc))
                    self.ax.plot(xSeq, meanCurveAux, linestyle='dotted')

            self.ax.plot(xSeq, meanCurve, label=label)

            if (vc > 0):
                for i in range(len(meanCurve)):
                    meanCurveAux[i] = meanCurve[i] * (1 - vc)
                self.ax.plot(xSeq, meanCurveAux, linestyle='dotted')

        if printLegend:
            self.ax.legend(loc="center left", bbox_to_anchor=[self.left, self.top],
                       ncol=self.legColumns, shadow=True, title=self.legendTitle)
            self.ax.get_legend().get_title().set_color("blue")

        
    def multiLinePrint(self):
        plt.draw()
        plt.show()       
                            
class Histogram_FreqDistribution:
    def __init__(self, desk, title):
        plt.close("all")
        plt.clf()
        
        if desk:
            self.fig = plt.figure(1, dpi=desk.dpi)
            self.desk = desk
            self.fontsize = 24- 2* int(round( (desk.dpi-100)/ 100.,0))
        else:
            self.fig = plt.figure(1, dpi=300) 
            self.desk = None
            self.fontsize = 15
        
        self.lines = 1
        self.columns = 2

        self.left=0.07
        self.right=0.95
        self.bottom=0.08
        self.top=0.92

        self.bar = None
        self.title = title
        self.myPlot = Plot()
        
        self.meanY = None
        self.medianY = None
        self.stdY  = None
        self.maxY = None
        self.minY = None
                            
                            
    def round2_20(self, val):
        if val > 200:
            return round(val,0)
        
        if val > 20:
            return round(val,1)
        
        if val > 2:
            return round(val,2)

        return round(val,4)

    
    def plot_H_MI(self, HorVert, xSeq, ySeq, ySE, showError=False, xminPar=None, yMinPar=None, yMaxPar=None):        
        if not xminPar:
            xMin = 0
        else:        
            xMin = xminPar

        if self.desk.mnat:
            unit = 'mnat'
        else:
            unit = 'nat'    
  
        xMax = len(ySeq)+1
        yMin = np.min([x for x in ySeq if x != 0])
        yMax = np.max([x for x in ySeq if x != 0])
        
        ySEmax = 0

        if showError and len(ySE) > 0:
            ySEmax = np.max(ySE)

        yMin -= ySEmax
        yMax += ySEmax
        
        if yMinPar:
            yMin = yMinPar
        else:
            if yMin < 0:
                yMin *= 1.05
            else:
                yMin *= 0.95
            
        if yMaxPar:
            yMax =yMaxPar
        else:
            if yMax + ySEmax < 0:
                yMax *= 0.95
            else:
                yMax *= 1.05 

        if HorVert == 'H':
            xTitle = 'k'
            if self.desk.isLog:
                yTitle = 'log(<MI>) (%s)'%(unit)
            else:
                yTitle = '<MI> (%s)'%(unit)
        else:
            xTitle = 'nuc'
            if self.desk.isLog:
                yTitle = 'log(<H>) (%s)'%(unit)
            else:
                yTitle = '<H> (%s)'%(unit)

        yLine = []

        if (showError):
            for _ in range(xMax):
                yLine.append(self.meanY)

        if self.desk.rand_method_var.get() == "shuffle" or self.desk.rand_method_var.get() == "random":
            roundVal = 3
        else:
            if self.desk.mnat:
                roundVal = 2
            else:
                roundVal = 1
                                                  
        titleG = self.title
        strfloat='.'+str(roundVal)
        tit = '= %xxxf(%xxxf) range = [%xxxf,%xxxf] %s' 
        tit = tit.replace('xxx', strfloat)                

        self.printSpectra(HorVert, xSeq, ySeq, ySE, yLine, xMin, xMax, yMin, yMax, xTitle, yTitle, titleG, 30, 'black')

    
    def printSpectra(self, HorVert, xSeq, ySeq, ySE, yLine, par_xmin, par_xmax, par_yMin, par_ymax, par_xlabel, par_ylabel, par_title, par_width, par_color):
        ax = plt.subplot(self.lines, self.columns, 1)
        plt.subplots_adjust(left=self.left, right=self.right, bottom=self.bottom, top=self.top)

        par_xmax = np.max(xSeq)
        plt.plot(xSeq, ySeq, 'r', color='green')
        plt.errorbar(x=xSeq, y=ySeq, yerr=ySE, ecolor='darkorange')

        stdY2 = 2 * self.stdY
                        
        seqSup = []
        seqMean = []
        
        if HorVert == 'H':
            seqInf = []
            for _ in range(len(xSeq)):
                seqSup.append(self.meanY+stdY2)
                seqMean.append(self.meanY)
                seqInf.append(self.meanY-stdY2)
            
            plt.plot(xSeq, seqSup, color='red')
            plt.plot(xSeq, seqMean, color='black')
            plt.plot(xSeq, seqInf, color='red')
        else:
            for _ in range(len(xSeq)):
                seqSup.append(self.meanY+stdY2)
                seqMean.append(self.meanY)
            
            plt.plot(xSeq, seqSup, color='red')
            plt.plot(xSeq, seqMean, color='black')


        ax.annotate(r'$2\sigma$', xy=(xSeq[10], self.meanY+stdY2*1.05), color='red')
                
        plt.xlim(par_xmin, par_xmax)
        plt.ylim(par_yMin, par_ymax)

        plt.title(par_title, fontsize=self.fontsize)
        plt.xlabel(par_xlabel, fontsize=self.fontsize)
        plt.ylabel(par_ylabel, fontsize=self.fontsize)
        plt.tick_params(labelsize=self.fontsize-2)
        

    def sameBar(self, xSeq, ySeq, linestyleCode, color):
        plt.subplot(1, 2, 1)

        yMax = np.max(ySeq)
        yMin = np.min(ySeq)
        
        if yMax > self.yMax:
            self.yMax = yMax
            
            if self.yMax >= 0:
                self.yMax *= 1.1
            else:
                self.yMax *= .9
                
        if yMin < self.yMin:
            self.yMin = yMin
            
            if self.yMin >= 0:
                self.yMin *= .9
            else:
                self.yMin *= 1.1
                            
        plt.ylim(self.yMin, self.yMax)        
        plt.plot(xSeq, ySeq, linestyleCode, color=color)
        

    def nextBar(self, xSeq, ySeq, yLine, parError, par_xmin, par_xmax, par_yMin, par_ymax, par_xlabel, par_ylabel, par_title, par_width, par_color):
        ax = plt.subplot(self.lines, self.columns, 2)
         
        wid = par_width
        
        if (len(yLine) == 0):
            # plt.bar(xSeq, ySeq, yerr=parError, width=wid, error_kw=dict(elinewidth=wid / 3, color=par_color, ecolor='red'))
            plt.bar(xSeq, ySeq, width=wid, color=par_color)
        else:
            plt.bar(xSeq, ySeq, yerr=parError, width=wid, error_kw=dict(elinewidth=wid / 3, color=par_color, ecolor='red'))  

            plt.plot(xSeq, yLine, color='red')
            seqSup = []
            seqInf = []
            for i in range(len(yLine)):
                seqSup.append(yLine[i] + parError)
                seqInf.append(yLine[i] - parError)
            plt.plot(xSeq, seqSup, color='black')
            plt.plot(xSeq, seqInf, color='black')

        # vertical line from (70,100) to (70, 250)
        if self.meanY:
            plt.plot([self.meanY, self.meanY], [par_yMin, par_ymax], 'k-', lw=2, color='black')
            plt.plot([self.meanY+self.stdY, self.meanY+self.stdY], [par_yMin, par_ymax], '--', lw=2, color='red')
            plt.plot([self.meanY-self.stdY, self.meanY-self.stdY], [par_yMin, par_ymax], '--', lw=2, color='red')
            plt.plot([self.medianY, self.medianY], [par_yMin, par_ymax], 'k-', lw=2, color='yellow')
            
            xPosAnn = self.meanY+self.stdY + .2*wid

            ax.annotate(r'$1\sigma$', xy=(xPosAnn, (par_ymax-par_yMin)*.75), color='red') #oi

        if ((par_yMin < 0.5) and (par_yMin > 0)):
            par_yMin = 0
            
        plt.xlim(par_xmin, par_xmax)
        plt.ylim(par_yMin, par_ymax)
        plt.title(par_title, fontsize=self.fontsize)
        
        plt.xlabel(par_xlabel, fontsize=self.fontsize)
        plt.ylabel(par_ylabel, fontsize=self.fontsize)
        plt.tick_params(labelsize=self.fontsize-2) 
                
                                
    def densityBar(self, HorVert, seq, bins=20, par_color="blue"):        
        plt.subplot(self.lines, self.columns, 2)
        
        if self.desk.mnat:
            unit = 'mnat'
        else:
            unit = 'nat'        

        ylabel = 'frequency (%)'
        if HorVert == 'H':
            if self.desk.isLog:
                xlabel = 'log(<MI>) (%s)'%(unit)
            else: 
                xlabel = '<MI> (%s)'%(unit)
        else:
            if self.desk.isLog:
                xlabel = 'log(<H>) (%s)'%(unit)
            else: 
                xlabel = '<H> (%s)'%(unit)
            
        meanY = np.mean(seq)
        stdY = np.sqrt(np.var(seq))
        medianY = np.median(seq)
        
        if self.desk.rand_method_var.get() == "shuffle" or self.desk.rand_method_var.get() == "random":
            roundVal = 3
        else:
            if self.desk.mnat:
                roundVal = 2
            else:
                roundVal = 1
                        
        strfloat='.'+str(roundVal)
        tit = 'Frequency Distribution \n mean=%xxxf(%xxxf); median=%xxxf %s'
        tit = tit.replace('xxx', strfloat)
        titleG = tit %(meanY, stdY, medianY, unit)
        
        
        ax = plt.subplot(self.lines, self.columns, 2)
        (n, _, _) = ax.hist(seq,bins=bins,color=par_color)
        # (n, bins, patches) = ax.hist(seq,bins=bins,color=par_color)

        yMax = np.max(n)
        yMax2 = yMax*1.05

        plt.plot([self.meanY, self.meanY], [0, yMax2], 'k-', lw=2, color='black')
        plt.plot([self.meanY+self.stdY, self.meanY+self.stdY], [0, yMax2], '--', lw=2, color='red')
        plt.plot([self.meanY-self.stdY, self.meanY-self.stdY], [0, yMax2], '--', lw=2, color='red')
        plt.plot([self.medianY, self.medianY], [0, yMax2], 'k-', lw=2, color='yellow')
        
        xPosAnn = self.meanY+self.stdY*1.05

        ax.annotate(r'$1\sigma$', xy=(xPosAnn, yMax*.75), color='red') #oi
            
        plt.ylim(0, yMax2)
        plt.title(titleG, fontsize=self.fontsize)
        
        plt.xlabel(xlabel, fontsize=self.fontsize)
        plt.ylabel(ylabel, fontsize=self.fontsize)
        plt.tick_params(labelsize=self.fontsize-2)
        
        return

        
    def printBar(self):
        plt.show()
        
    def savePlot(self, desk, pictureName):
        plt.savefig(pictureName, format=desk.imgType, dpi=desk.dpi)
        
class HeatMap:
    # tirei o save e criei o print
    def __init__(self, desk, is3D, lenSeq, tickWidth=50):
        self.desk = desk
        
        ''' http://stackoverflow.com/questions/22408237/named-colors-in-matplotlib '''
        self.colorList = ['darkblue', 'lightblue', 'g', 'purple', 'r']
        self.markerList = ['.','1','2','D','o']
        
        plt.close("all")
        plt.clf()
        
        if is3D:
            # Twice as wide as it is tall.
            # self.fig = plt.figure(figsize=plt.figaspect(0.5), dpi=desk.dpi)
            self.fig = plt.figure(1, dpi=desk.dpi)
        else:
            self.fig = plt.figure(1, figsize=(8, 8), dpi=desk.dpi)
        
        self.lines = 1
        self.columns = 1

        self.fontsize = 24- 2* int(round( (desk.dpi-100)/ 100.,0))
        
        self.left = 0.10
        self.bottom = 0.10
        self.width = 0.80
        self.height = .75

        
        self.dna_prot = desk.dna_prot
        self.dpi = desk.dpi
        
        self.tick = [i*tickWidth for i in range(lenSeq / tickWidth)]
  
        self.myPlot = Plot()
            
          
    def colors(self):
        pass
        '''
        cnames = {
        'aliceblue':            '#F0F8FF',
        'antiquewhite':         '#FAEBD7',
        'aqua':                 '#00FFFF',
        'aquamarine':           '#7FFFD4',
        'azure':                '#F0FFFF',
        'beige':                '#F5F5DC',
        'bisque':               '#FFE4C4',
        'black':                '#000000',
        'blanchedalmond':       '#FFEBCD',
        'blue':                 '#0000FF',
        'blueviolet':           '#8A2BE2',
        'brown':                '#A52A2A',
        'burlywood':            '#DEB887',
        'cadetblue':            '#5F9EA0',
        'chartreuse':           '#7FFF00',
        'chocolate':            '#D2691E',
        'coral':                '#FF7F50',
        'cornflowerblue':       '#6495ED',
        'cornsilk':             '#FFF8DC',
        'crimson':              '#DC143C',
        'cyan':                 '#00FFFF',
        'darkblue':             '#00008B',
        'darkcyan':             '#008B8B',
        'darkgoldenrod':        '#B8860B',
        'darkgray':             '#A9A9A9',
        'darkgreen':            '#006400',
        'darkkhaki':            '#BDB76B',
        'darkmagenta':          '#8B008B',
        'darkolivegreen':       '#556B2F',
        'darkorange':           '#FF8C00',
        'darkorchid':           '#9932CC',
        'darkred':              '#8B0000',
        'darksalmon':           '#E9967A',
        'darkseagreen':         '#8FBC8F',
        'darkslateblue':        '#483D8B',
        'darkslategray':        '#2F4F4F',
        'darkturquoise':        '#00CED1',
        'darkviolet':           '#9400D3',
        'deeppink':             '#FF1493',
        'deepskyblue':          '#00BFFF',
        'dimgray':              '#696969',
        'dodgerblue':           '#1E90FF',
        'firebrick':            '#B22222',
        'floralwhite':          '#FFFAF0',
        'forestgreen':          '#228B22',
        'fuchsia':              '#FF00FF',
        'gainsboro':            '#DCDCDC',
        'ghostwhite':           '#F8F8FF',
        'gold':                 '#FFD700',
        'goldenrod':            '#DAA520',
        'gray':                 '#808080',
        'green':                '#008000',
        'greenyellow':          '#ADFF2F',
        'honeydew':             '#F0FFF0',
        'hotpink':              '#FF69B4',
        'indianred':            '#CD5C5C',
        'indigo':               '#4B0082',
        'ivory':                '#FFFFF0',
        'khaki':                '#F0E68C',
        'lavender':             '#E6E6FA',
        'lavenderblush':        '#FFF0F5',
        'lawngreen':            '#7CFC00',
        'lemonchiffon':         '#FFFACD',
        'lightblue':            '#ADD8E6',
        'lightcoral':           '#F08080',
        'lightcyan':            '#E0FFFF',
        'lightgoldenrodyellow': '#FAFAD2',
        'lightgreen':           '#90EE90',
        'lightgray':            '#D3D3D3',
        'lightpink':            '#FFB6C1',
        'lightsalmon':          '#FFA07A',
        'lightseagreen':        '#20B2AA',
        'lightskyblue':         '#87CEFA',
        'lightslategray':       '#778899',
        'lightsteelblue':       '#B0C4DE',
        'lightyellow':          '#FFFFE0',
        'lime':                 '#00FF00',
        'limegreen':            '#32CD32',
        'linen':                '#FAF0E6',
        'magenta':              '#FF00FF',
        'maroon':               '#800000',
        'mediumaquamarine':     '#66CDAA',
        'mediumblue':           '#0000CD',
        'mediumorchid':         '#BA55D3',
        'mediumpurple':         '#9370DB',
        'mediumseagreen':       '#3CB371',
        'mediumslateblue':      '#7B68EE',
        'mediumspringgreen':    '#00FA9A',
        'mediumturquoise':      '#48D1CC',
        'mediumvioletred':      '#C71585',
        'midnightblue':         '#191970',
        'mintcream':            '#F5FFFA',
        'mistyrose':            '#FFE4E1',
        'moccasin':             '#FFE4B5',
        'navajowhite':          '#FFDEAD',
        'navy':                 '#000080',
        'oldlace':              '#FDF5E6',
        'olive':                '#808000',
        'olivedrab':            '#6B8E23',
        'orange':               '#FFA500',
        'orangered':            '#FF4500',
        'orchid':               '#DA70D6',
        'palegoldenrod':        '#EEE8AA',
        'palegreen':            '#98FB98',
        'paleturquoise':        '#AFEEEE',
        'palevioletred':        '#DB7093',
        'papayawhip':           '#FFEFD5',
        'peachpuff':            '#FFDAB9',
        'peru':                 '#CD853F',
        'pink':                 '#FFC0CB',
        'plum':                 '#DDA0DD',
        'powderblue':           '#B0E0E6',
        'purple':               '#800080',
        'red':                  '#FF0000',
        'rosybrown':            '#BC8F8F',
        'royalblue':            '#4169E1',
        'saddlebrown':          '#8B4513',
        'salmon':               '#FA8072',
        'sandybrown':           '#FAA460',
        'seagreen':             '#2E8B57',
        'seashell':             '#FFF5EE',
        'sienna':               '#A0522D',
        'silver':               '#C0C0C0',
        'skyblue':              '#87CEEB',
        'slateblue':            '#6A5ACD',
        'slategray':            '#708090',
        'snow':                 '#FFFAFA',
        'springgreen':          '#00FF7F',
        'steelblue':            '#4682B4',
        'tan':                  '#D2B48C',
        'teal':                 '#008080',
        'thistle':              '#D8BFD8',
        'tomato':               '#FF6347',
        'turquoise':            '#40E0D0',
        'violet':               '#EE82EE',
        'wheat':                '#F5DEB3',
        'white':                '#FFFFFF',
        'whitesmoke':           '#F5F5F5',
        'yellow':               '#FFFF00',
        'yellowgreen':          '#9ACD32'}
        '''
        
        
    def plotMI(self, title, seq, limSup):
        bottom = 0.08
        
        ax1 = plt.subplot(1, 2, 1)
        plt.subplots_adjust(left=0.05, right=0.97, bottom=bottom, top=0.82, wspace=0.3)

        if limSup:
            # , cmap=plt.cm.YlGnBu   cmap=plt.cm.get_cmap('OrRd'),
            im = ax1.matshow(seq, aspect='auto', origin='lower', cmap=plt.cm.get_cmap(self.desk.heatmap_color), vmin=0, vmax=limSup) 
        else:
            im = ax1.matshow(seq, aspect='auto', origin='lower') # , cmap=plt.cm.YlGnBu

        ax1.set_xticks(self.tick)
        ax1.set_yticks(self.tick)

        #plt.title("")
        self.fig.text(.5, .92, title, ha='center', fontsize=self.fontsize, color="black")

        ax1.set_xticklabels(self.tick,minor=False,fontsize=self.fontsize-4)
        ax1.set_yticklabels(self.tick,minor=False,fontsize=self.fontsize-4)        
        
        if self.dna_prot=='DNA':
            label = 'bp-nuc'
            plt.xlabel(label, fontsize=self.fontsize)
            plt.ylabel(label, fontsize=self.fontsize)
        else:
            label = 'aa'
            plt.xlabel(label, fontsize=self.fontsize)
            plt.ylabel(label, fontsize=self.fontsize)
            
        
        width = 0.02
        height =.74

        axcolor = self.fig.add_axes([.46, bottom, width, height])
        cbar = plt.colorbar(im, cax=axcolor)
        cbar.ax.tick_params(labelsize=self.fontsize-2)
        plt.axhline()
    

    def densityHeatmapBar(self, seq, limSup, bins=20, par_color="blue"):
        if self.desk.mnat:
            unit = 'mnat'
        else:
            unit = 'nat'        

        ylabel = 'frequency (log10(#)'
        
        if self.desk.isLog:
            xlabel = 'log(<VMI>) (%s)'%(unit)
        else: 
            xlabel = '<H> (%s)'%(unit)
        
        meanY = np.mean(seq)
        stdY = np.sqrt(np.var(seq))
        medianY = np.median(seq)
      
        if self.desk.rand_method_var.get() == "shuffle" or self.desk.rand_method_var.get() == "random":
            roundVal = 3
        else:
            if self.desk.mnat:
                roundVal = 2
            else:
                roundVal = 1
                        
        strfloat='.'+str(roundVal)
        tit = 'Frequency Distribution \n mean=%xxxf(%xxxf); median=%xxxf %s'
        tit = tit.replace('xxx', strfloat)
        titleG = tit %(meanY, stdY, medianY, unit)
        

        ax = plt.subplot(1, 2, 2)
        (n, _, _) = ax.hist(seq,bins=bins,color=par_color,log=True)
        # (n, bins, patches) = ax.hist(seq,bins=bins,color=par_color)

        yMax = np.max(n)
        yMax2 = yMax*1.05

        
        plt.plot([meanY, meanY], [0, yMax2], 'k-', lw=2, color='black')
        plt.plot([meanY+stdY, meanY+stdY], [0, yMax2], '--', lw=2, color='red')
        plt.plot([meanY-stdY, meanY-stdY], [0, yMax2], '--', lw=2, color='red')
        plt.plot([medianY, medianY], [0, yMax2], 'k-', lw=2, color='yellow')
        
        xPosAnn = meanY+stdY*1.05

        ax.annotate(r'$1\sigma$', xy=(xPosAnn, yMax*.75), color='red') #oi
            
        plt.ylim(0, yMax2)
        plt.xlim(0, limSup)
        
        plt.title(titleG, fontsize=self.fontsize)
        
        plt.xlabel(xlabel, fontsize=self.fontsize)
        plt.ylabel(ylabel, fontsize=self.fontsize)
        plt.tick_params(labelsize=self.fontsize-2)
        
        return 

    def plotMI_3D(self, title, seqX, seqY, seqZ, L, limSup, unit):

        #---- First subplot
        ax = self.fig.add_subplot(1, 1, 1, projection='3d')

        #ax.set_xticklabels(self.tick,fontsize=12)
        #ax.set_yticklabels(self.tick,fontsize=12) 

        plt.xticks(np.arange(0, L+1, 50), rotation='vertical')
        plt.yticks(np.arange(0, L+1, 50), rotation='vertical')
        '''
        print 'l',L
        L = 1600
        aX=[0,L]
        aY=[0,L]
        aZ=[20,300]
        ax.scatter(aX, aY, aZ, c='red', marker='o')
        '''
        for i in range(5):
            ax.scatter(seqX[i], seqY[i], seqZ[i], c=self.colorList[i], marker=self.markerList[i])

        '''
        aX=[20,L-20]
        aY=[20,L-20]
        aZ=[50,250]        
        ax.scatter(aX, aY, aZ, c='blue', marker='o')
        '''
        ax.set_xlabel('nuc', fontsize=12)
        ax.set_ylabel('nuc', fontsize=12)
        ax.set_zlabel('<VMI> (' + unit + ')', fontsize=12)

        ax.set_xlim3d(0, L)
        ax.set_ylim3d(0, L)
        ax.set_zlim3d(0, limSup)

        plt.title(title, fontsize=14)


    def printBar(self):
        plt.show()
        
    def savePlot(self, desk, pictureName):
        plt.savefig(pictureName, format=desk.imgType, dpi=desk.dpi)

        
class PrintDendogram:
    def __init__(self, seq, names):
       
        """
        linkage(y, method='single', metric='euclidean'):
    
        Performs hierarchical/agglomerative clustering on the
        condensed distance ma  number of original observations paired
        in the distance matrix. The behavior of this function is very
        similar to the MATLAB(TM) linkage function. 
        """
        
        Z = linkage(seq)
        
        d = dendrogram(Z)
        
        d['leaves'] = names
        plt.show()

        
               
class HistogramGraphic:
    def __init__(self, fontsize=10):
        self.fontsize = fontsize
        
    def gHist(self, numFig, num, seq, par_xlabel, par_ylabel, par_title, par_color="blue", bins=20):
        plt.subplot(1, numFig, num)
        plt.hist(seq, bins=bins, color=par_color)
        '''
        (n, bins, patches) = plt.hist(parY)
        
        print 'n', n
        print 'bins', bins
        print 'patches', patches
        '''
        plt.ylabel(par_ylabel, fontsize=self.fontsize-2)
        plt.xlabel(par_xlabel, fontsize=self.fontsize-2)
        plt.title(par_title, fontsize=self.fontsize)
        
    def printBar(self):
        plt.show()
        
                
class PlotGraphic:
    def __init__(self, num, par_xseq, par_ySeq, par_xmin, par_xmax, par_yMin, par_ymax, par_xlabel, par_ylabel, par_title):
        xseq = []
        ySeq = []
        xMin = []
        xMax = []
        yMin = []
        yMax = []
        xlabel = []
        ylabel = []
        title = []
        figure = []

        self.xseq = xseq
        self.ySeq = ySeq
        self.xMin = xMin
        self.xMax = xMax
        self.yMin = yMin
        self.yMax = yMax
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.title = title
        
        self.figure = figure

        self.init(num, par_xseq, par_ySeq, par_xmin, par_xmax, par_yMin, par_ymax, par_xlabel, par_ylabel, par_title)
        
        # print "inicializing ...", num, " '" + self.title[num] + "'"  
        
    def init(self, num, par_xseq, par_ySeq, par_xmin, par_xmax, par_yMin, par_ymax, par_xlabel, par_ylabel, par_title):
        self.xseq.append(par_xseq)
        self.ySeq.append(par_ySeq)
        self.xMin.append(par_xmin)
        self.xMax.append(par_xmax)
        self.yMin.append(par_yMin)
        self.yMax.append(par_ymax)
        self.xlabel.append(par_xlabel)
        self.ylabel.append(par_ylabel)
        self.title.append(par_title)
        
        # self.figure.append(plt.figure(num))
        # plt.subplot(1, 2, num)

        # print "inicializing ...", num, " '" + self.title[num] + "'"  
        
        
        
    def buildPlot(self, num, lineColor, backColor):
        fig = self.figure[num]
        
        ax = fig.add_subplot(111)
        # the histogram of the data
        ax.hist(self.xseq[num], color=lineColor, facecolor=backColor) 

        
        # l = ax.plot(bincenters, y, 'r--', linewidth=1)
        ax.plot(self.xseq[num], self.ySeq[num], 'r--', linewidth=1)
        
        ax.set_xlabel(self.xlabel[num])
        ax.set_ylabel(self.ylabel[num])
        ax.set_title(self.title[num])
        
        # print 'x length ', len(self.xseq[num])
        ax.set_xlim(self.xMin[num], self.xMax[num])
        # print 'np.max value ',  np.max(self.ySeq)
        ax.set_ylim(self.yMin[num], self.yMax[num])

        ax.grid(True)
        

    
    def printPlot(self):
        plt.show()        


class GeneralGraphic:
    def __init__(self, desk):
        self.myPlot = Plot()
        self.fontsize = 8
        self.fig = plt.figure(1, dpi=desk.dpi)
        plt.subplots_adjust(left=0.1, right=0.95, bottom=0.15, top=0.90)
        
        # print_graph(self, desk, fig, pictureName, frame=None, stay=False):

    def histogram_H0_Ha(self, desk, listH0, listHa, label_random):
        maxi = np.max(listH0)
        maxi2 = np.max(listHa)
            
        if maxi2 > maxi:
            maxi = maxi2

        # fig = plt.figure(1, dpi=desk.dpi)
        ax = self.fig.add_subplot("122")

        ax.hist(listH0, bins=20, color="blue",alpha=0.3, label='H0: same species', edgecolor = "blue")
        ax.hist(listHa, bins=20, color="red", alpha=0.5, label='Ha: diff. species', edgecolor = "red")
        
        plt.ylabel("frequency", fontsize=self.fontsize-2)
        plt.xlabel("HMI in %s"%(desk.unit), fontsize=self.fontsize-2)
        plt.xlim(0, maxi)
        plt.legend(loc='upper right', fontsize=self.fontsize-5)
        if label_random != "":
            label_random = " " + label_random
        plt.title("JSD[HMI(%s%s)]%s"%(desk.minmax,desk.str_correction,label_random), fontsize=self.fontsize-1)
        plt.tick_params(labelsize=self.fontsize-2)

        
    def ROC_Curve(self, desk, sensitivityList, specificityList, label_random):
        y = np.array(sensitivityList)
        x = 1.-np.array(specificityList)

        ax = self.fig.add_subplot("121")

        ax.scatter(x,y, edgecolor = "none")
        plt.ylabel("TPR = sensitivity", fontsize=self.fontsize-2)
        plt.xlabel("FPR = 1-specificity", fontsize=self.fontsize-2)
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        
        if label_random != "":
            label_random = " " + label_random
       
        plt.title("ROC for HMI(%s%s)%s"%(desk.minmax,desk.str_correction,label_random), fontsize=self.fontsize-1)
        
        ax.tick_params(labelsize=self.fontsize-2)
        

        
        
        