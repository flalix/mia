#-*- coding: utf-8 -*-
'''
Created on Jul 29, 2015
Updated on Oct 10, 2015

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica
'''
import classes.mrBayesClass as mb
import numpy as np
import matplotlib.pyplot as plt
import classes.BarGraphic as graphPac

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

        mbs = mb.mrBayesClass(desk)
        
        if not mbs.read_runPs():
            self.error_msg = 'Problems reading Mr.Bayes files'
            return

        if not mbs.stat_files_to_dic():
            self.error_msg = 'Could not read Mr.Bayes statistic files'
            return
        
        mbs.summary = {}
        
        params = []
        if desk.piA_var.get(): params.append("pi(A)")
        if desk.piC_var.get(): params.append("pi(C)")
        if desk.piG_var.get(): params.append("pi(G)")
        if desk.piT_var.get(): params.append("pi(T)")
        if desk.rAC_var.get(): params.append("r(A<->C)")
        if desk.rAG_var.get(): params.append("r(A<->G)")
        if desk.rAT_var.get(): params.append("r(A<->T)")
        if desk.rCG_var.get(): params.append("r(C<->G)")
        if desk.rCT_var.get(): params.append("r(C<->T)")
        if desk.rGT_var.get(): params.append("r(G<->T")
        if desk.LnL_var.get(): params.append("LnL")
        if desk.LnPr_var.get(): params.append("LnPr")
        if desk.TL_var.get(): params.append("TL")
        if desk.alpha_var.get(): params.append("alpha")
        if desk.off_on_var.get(): params.append("s(off->on)")
        if desk.on_off_var.get(): params.append("s(on->off)")
        if desk.pinvar_var.get(): params.append("pinvar")

        if not params:  
            self.error_msg = 'Define at least one param.'
            return
        
        iPar = 0
        
        numCols = 46
        numLines = 16
        self.myPlot = graphPac.Plot()
        
        for par in params:
            iPar += 1
            mbs.summary[par] = {}
            
            #fig = figList[iPar]
            fig = plt.figure(iPar)
            plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.90)
               
            # print plt.get_backend()
            mng = plt.get_current_fig_manager()
            mng.window.wm_geometry("1400x900+50+50")

            seqXs = []
            mini, maxi = float('inf'), -float('inf')
            
            for min_max in mbs.min_maxs:
                for alig_cons in mbs.alig_conss:
                    for std_covar in mbs.std_covars: 
                           
                        study = min_max + '-' + alig_cons + '-' + std_covar
                     
                        try:
                            vals = np.array(mbs.dic_detaild_params[study][par])

                            if np.min(vals) < mini:
                                mini = min(vals)
                            elif np.max(vals) > maxi:
                                maxi = max(vals)                  
                        except:
                            pass
                        
            if par != "LnL" and par != "LnPr":
                if par == "alpha":
                    mini = 0
                else:
                    if mini < 0: mini = 0
            
            if par in ["pi(A)","pi(C)","pi(G)","pi(T)"]:
                if maxi > 1: maxi = 1

            ticks = []
            ticks0 = []
            div = 5
            delta = (maxi-mini)/float(div-1)
            for i in range(div):
                ticks.append(round(mini+delta*i,2))
                ticks0.append(0)
                                    
            
            cntCols = 0
            stri_summary = "   \tAligned\tConsensus\n"
            for min_max in mbs.min_maxs:
                for alig_cons in mbs.alig_conss:
                    
                    listData = []
                   
                    cntLine = 0
                    
                    for std_covar in mbs.std_covars: 
                        try:
                            _ = mbs.mb_filenames[min_max][alig_cons][std_covar]
                            didFind = True
                        except:
                            didFind = False
                            continue
                                             
                           
                        study = min_max + '-' + alig_cons + '-' + std_covar
                        
                        '''
                        try:
                            dic = mbs.dic_detaild_params[study]
                        except:
                            # if file not found or could not read: it is not in dic_detailed_params 
                            continue
                        '''
                        
                        mbs.summary[par][study] = {}
                      
                        try:
                            vals = np.array(mbs.dic_detaild_params[study][par])
                            
                            mbs.summary[par][study] = {}
                            N = len(vals)
                            seqXs.append(np.array(vals))
        
                           
                            mu, sigma, hmu, hsigma, q025, _, q2, _, q975 = mbs.stat_ppf_decimal(vals, par)

                            mbs.summary[par][study]["N"] = N
                            mbs.summary[par][study]["mu"] = mu
                            mbs.summary[par][study]["sigma"] = sigma
                            SE = sigma/np.sqrt(N)
                            mbs.summary[par][study]["SE"] = SE
                            mbs.summary[par][study]["N"] = N
                            mbs.summary[par][study]["median"] = q2
                            mbs.summary[par][study]["q025"] = q025
                            mbs.summary[par][study]["q975"] = q975
                            
                            if hmu:
                                mbs.summary[par][study]["hmu"] = hmu
                                mbs.summary[par][study]["hsigma"] = hsigma
                                hSE = hsigma/np.sqrt(N)
                                mbs.summary[par][study]["hSE"] = hSE
                                
                            mbs.summary[par][study]["min"] = np.min(vals)
                            mbs.summary[par][study]["max"] = np.max(vals)


                            try:
                                ''' LnL doesn't have ESS and PSRF - adopted from TL
                                    total tree length (the sum of all branch lengths, TL) '''
                                if par == "LnL":
                                    mbs.summary[par][study]["avgESS"] = mbs.dic_pstat[study]["TL"]['avgESS']
                                else:
                                    mbs.summary[par][study]["avgESS"] = mbs.dic_pstat[study][par]['avgESS']  
                            except:
                                mbs.summary[par][study]["avgESS"] = 0
                                
                            try:
                                if par == "LnL":
                                    mbs.summary[par][study]["PSRF"] = mbs.dic_pstat[study]["TL"]['PSRF']
                                else:
                                    mbs.summary[par][study]["PSRF"] = mbs.dic_pstat[study][par]['PSRF']
                            except:
                                mbs.summary[par][study]["PSRF"] = 10000
                        
                            if mbs.summary[par][study]["avgESS"] < 90:
                                sAvgESS = "avgESS=%3.1f **"%mbs.summary[par][study]["avgESS"]
                            else:
                                sAvgESS = "avgESS=%3.1f"%mbs.summary[par][study]["avgESS"]
                                
                            if (mbs.summary[par][study]["PSRF"] - 1) > .1:
                                sPSRF = " PSRF=%1.2f ***"%mbs.summary[par][study]["PSRF"]
                            else:
                                sPSRF = " PSRF=%1.2f"%mbs.summary[par][study]["PSRF"]
        
        
                            if cntLine == 0:
                                ax = plt.subplot2grid((numLines,numCols), (0, cntCols), rowspan=4,colspan=10)
                            else:
                                ax = plt.subplot2grid((numLines,numCols), (6, cntCols), rowspan=4,colspan=10)
        
                            if par == "LnL":
                                print par, study, mu, sigma, hmu, hsigma
                                stri_summary = mbs.show_lnl(desk, plt, ax, study, sAvgESS, sPSRF, cntLine, cntCols, min_max, alig_cons, stri_summary)
                            else:
                                print par, study, mu, sigma
                                mbs.show_distrib(plt, ax, par, study, sAvgESS, sPSRF, ticks, desk.colors, iPar, mini, maxi, cntCols)
                          
                            listData.append(vals) 
                            cntLine += 1
                        except:
                            pass
        
                    if didFind:    
                        if len(listData) == 2:
                            if par == "LnL":  
                                mbs.show_LRT(alig_cons, par, plt, numLines,numCols, cntCols)
                            else:
                                mbs.show_conf_interval(min_max, alig_cons, par, plt, numLines,numCols, cntCols, ticks, ticks0, mini, maxi, listData)
                            
                            mbs.show_qqplot2(min_max, alig_cons, par,plt, numLines,numCols, cntCols, listData)

                        else:
                            pass

                        cntCols += 12 
            
            if par == "LnL":
                stri = ""
            else:
                f_value, p_value = mbs.calc_anova(seqXs)
                if p_value <= 0.05:
                    stri = ', at least one distribution is statistically different.'
                else:
                    stri =  ', the distributions are statistically similar.'
                    
                stri += 'ANOVA: f-value %2.3e   p_value %2.3e' %(f_value, p_value)

            left = .05
            top = .97
            fig.text(left, top, par+stri, color="red") 
            
            stri = mbs.ttest_summary(par)
            print stri
            
            sPar = par.replace(">","").replace("<","")
            pictureName = "%s_%s_mr_bayes_analisis_param_%s"%(desk.organism, desk.gene_title, sPar)
            self.myPlot.print_graph(self.desk, fig, pictureName, frame=self.desk.tk_root, stay=False)

        print stri_summary
        
        if desk.saveData:
            mbs.save_params(params)


        self.failure = False  
        return

                                            
          
                
                    