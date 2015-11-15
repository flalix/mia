#-*- coding: utf-8 -*-
'''
Created on Aug 27, 2015
Updated on Oct 10, 2015

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica
'''
# import FileDialog
import ete2 # Tree, PhyloTree, TreeStyle
import Bio.Phylo as Phylo
import os
import decimal
from copy import deepcopy
import numpy as np #, scipy.stats
from scipy.stats import f_oneway, chi2, norm, ttest_ind, probplot, hmean
from scipy.stats.mstats import mquantiles
import matplotlib

class mrBayesClass(object):
    '''
    MrBayes class:
        calc all MrBayes PD
        calc LRT
        for tree plotting:
            ete2 otherwise Phylo
            call: FigTree (external call) - build_fig_tree()
    '''
       
    def __init__(self, desk):
        self.desk = desk
        
        desk.isLog = False
        
        del_substring = '_%s'%(desk.organism)
        del_to_right = '_%s_Prod_'%(desk.gene_title)
        ''' which one: '= [&U] '   '= [&R] ' '''
        mb_clock = '= [&U] '

        self.del_substring = del_substring.lower()
        self.del_to_right = del_to_right.lower()
        self.mb_clock = mb_clock
        str_start = '   tree gen.'
        self.len_start = len(str_start)


        self.min_maxs = ['mincut', 'maxmer']
        self.alig_conss = ['aligned', 'consensus']
        self.std_covars = ['standard', 'covarion']

       
        self.dic_params = {}
        self.dic_pstat = {}
        self.dic_ll = {}
        
        ''' header run.tree has 5 lines '''
        self.num_lines_header_tree = 5
        ''' header run.p has 2 lines '''
        self.num_lines_header_pi = 2

        ''' dictionary result '''
        self.dic_generation = {}
        
        ''' list of Mr.Bayes files found '''
        self.mb_filenames = {}

        
    def read_runPs(self):
        self.dic_detaild_params = {}
        # self.dic_lnl = {}

        found = True
        for alig_cons in self.alig_conss:
            for min_max in self.min_maxs:
                for std_covar in self.std_covars:
                    
                    if std_covar == "standard":
                        if min_max == "mincut":
                            if alig_cons == "aligned":
                                filename = self.desk.mrBayes_fn01
                            else:
                                filename = self.desk.mrBayes_fn02
                        else:
                            if alig_cons == "aligned":
                                filename = self.desk.mrBayes_fn03
                            else:
                                filename = self.desk.mrBayes_fn04
                    else:
                        if min_max == "mincut":
                            if alig_cons == "aligned":
                                filename = self.desk.mrBayes_fn01_cov
                            else:
                                filename = self.desk.mrBayes_fn02_cov
                        else:
                            if alig_cons == "aligned":
                                filename = self.desk.mrBayes_fn03_cov
                            else:
                                filename = self.desk.mrBayes_fn04_cov
                    
                    filename = self.desk.mrBayes_path + filename
        
                    if os.path.isfile(filename):
                        if min_max not in self.mb_filenames.keys():
                            self.mb_filenames[min_max] = {}
                        dic = self.mb_filenames[min_max]
                            
                        if alig_cons not in dic.keys():
                            dic[alig_cons] = {}
                        dic2 = dic[alig_cons]
                            
                        dic2[std_covar] = filename
                        
                        if not self.read_run1p(filename, min_max, alig_cons, std_covar):
                            print('### Problems reading - %s'%(filename))
                        else:
                            print('>>> Ok - %s'%(filename))
                    else:
                        print('!!! Not Found - %s'%(filename))

                    
        return found

    def stat_ppf(self, x, par):

        x1, x2, x3, x4, x5 = mquantiles(x, prob=[0.025, 0.25, 0.5, 0.75, .975])
        mu = np.mean(x)
        sigma = np.sqrt(np.var(x))

        return mu, sigma, x1, x2, x3, x4, x5
                     
    def stat_ppf_decimal(self, x, par):

        if par == "LnL":
            # ini = int(len(x)*self.desk.burnin)+3
            # x = [decimal.Decimal(a).exp() for a in x[ini:]]
            # x1, x2, x3, x4, x5 = self.my_quantiles(x)
            x1, x2, x3, x4, x5 = mquantiles(x, prob=[0.025, 0.25, 0.5, 0.75, .975])
            mu = np.mean(x)
            sigma = np.sqrt(np.var(x))

            ini = int(len(x)*self.desk.burnin)+3
            x = [decimal.Decimal(a).exp() for a in x[ini:]]
            hmu, hsigma = self.harmonic_mean_sigma_decimal(x)
            
            return mu, sigma, round(hmu.ln(),1), -round(hsigma.ln(),1), x1, x2, x3, x4, x5
            
        else:
            x1, x2, x3, x4, x5 = mquantiles(x, prob=[0.025, 0.25, 0.5, 0.75, .975])
            mu = np.mean(x)
            sigma = np.sqrt(np.var(x))
        
        if par == "LnL":
            return mu, sigma, -hmu, hsigma, x1, x2, x3, x4, x5
            '''
            return round(mu.ln(),1), -round(sigma.ln(),1), round(hmu.ln(),1), -round(hsigma.ln(),1), \
                   round(x1.ln(),1), round(x2.ln(),1), round(x3.ln(),1), round(x4.ln(),1), round(x5.ln(),1)
            '''
        else:
            return mu, sigma, None, None, x1, x2, x3, x4, x5

    def my_quantiles(self, x):
        x = sorted(x)
        L = len(x)

        qls = []
        for alpha in [.025, .25, .5, .75, .975]:
            La = alpha*(L+1)
            diff = decimal.Decimal(La - np.floor(La))
            
            if int(np.floor(La)) < 1:
                qls.append( x[0])
            else:
                floor = int(np.floor(La)) - 1
                
                if floor >= L-1:
                    qls.append(x[L-1])
                else:
                    ceil = int(np.ceil(La)) - 1
                    qls.append(x[floor] + (x[ceil] - x[floor]) * diff)

        return qls[0], qls[1], qls[2], qls[3], qls[4] 
    
    def harmonic_mean_sigma(self, x):
        L = len(x)
        x = [1. / a for a in x]
        
        ''' harmonic mean '''
        mu = 1./np.mean(x)
        sigma = np.sqrt((np.mean(x))**(-4)*np.var(x)/L)
        
        return mu, sigma
    
        
    def harmonic_mean_sigma_decimal(self, x):
        L = len(x)
        x = [decimal.Decimal(1) / a for a in x]
        
        ''' harmonic mean '''
        mu = decimal.Decimal(1)/np.mean(x)
        ''' correct '''
        sigma = np.sqrt((np.mean(x))**(-4)*np.var(x)/L)
        
        return mu, sigma
    
    def my_harmonic_ppf(self, x):
        x = np.array(x)
        x = 1./x
        L = len(x)
        
        ''' harmonic mean '''
        mu = 1./np.mean(x)
        sigma = np.sqrt((np.mean(x))**(-4)*np.var(x)/L)

        x = sorted(x)
        qls = []
        for alpha in [.025, .25, .5, .75, .975]:
            La = alpha*(L+1)
           
            if int(np.floor(La)) < 1:
                qls.append( x[0])
            else:
                i = int( round(La))
                if i > L-1:
                    qls.append(x[L-1])
                else:
                    qls.append(x[i-1])

        return mu, sigma, 1./qls[4], 1./qls[3], 1./qls[2], 1./qls[1], 1./qls[0]
            
    def read_run1p(self, filename, min_max, alig_cons, std_covar):
        filename += '.run1.p'
        try:
            f = open(filename)
        except:
            print 'could not read %s'%(filename)
            return False
        
        study = min_max + '-' + alig_cons + '-' + std_covar
        self.dic_detaild_params[study] = {}
        
        lines = f.readlines()
        
        header = lines[1].rstrip().split('\t')

        ini = int(len(lines)*self.desk.burnin)+3
        
        for par in header:
            self.dic_detaild_params[study][par] = []
        
        for line in lines[ini:]:
            mat = line.rstrip().split('\t')
            if len(mat) != len(header):
                continue
            
            for i in range(len(mat)):
                par = header[i]
                try:
                    ''' generations: Gen '''
                    if i==0:
                        self.dic_detaild_params[study][par].append(int(mat[i]))
                    else:
                        ''' else: LnL    LnPr    TL    r(A<->C) .... '''
                        self.dic_detaild_params[study][par].append(float(mat[i]))
                except:
                    if i==0:
                        print "found some error converting Generation like integer",str(mat[i])
                    else:
                        print "found some error converting values like float",str(mat[i])
                
        '''
        self.dic_lnl[study] = {}
        self.dic_lnl[study]["Gen"] = []
        self.dic_lnl[study]["LnL"] = []
        
        if self.desk.isLog:
            for line in lines[ini:]:  #lines[10:]:
                mat = line.rstrip().split('\t')
                if len(mat) != len(header):
                    continue
    
                try:
                    self.dic_lnl[study]["Gen"].append(np.log10(int(mat[0])))
                except:
                    print "found some error in %s converting Generation like integer %s"%(study, str(mat[0]))
                    return False
                try:
                    if float(mat[1]) > 0:
                        print "found some error in %s converting LnL not negative %s"%(study, str(mat[1]))
                        return False
    
                    self.dic_lnl[study]["LnL"].append(np.log10(-float(mat[1])))
                except:
                    print "found some error converting LnL like float",str(mat[1])
        else:
            for line in lines[ini:]:  #lines[10:]:
                mat = line.rstrip().split('\t')
                if len(mat) != len(header):
                    continue
    
                try:
                    self.dic_lnl[study]["Gen"].append(np.log10(int(mat[0])))
                except:
                    print "found some error in %s converting Generation like integer %s"%(study, str(mat[0]))
                    return False
                try:
                    self.dic_lnl[study]["LnL"].append(float(mat[1]))
                except:
                    print "found some error converting LnL like float",str(mat[1])
        '''
        return True
    
            
    def calc_anova(self, seqXs):
        #one way anova
        f_value, p_value = f_oneway(*seqXs)
       
        return f_value, p_value
    
    
    def read_pstat(self, filename, min_max, alig_cons, std_covar):
        filename += '.pstat'
        try:
            f = open(filename)
        except:
            print 'could not read %s'%(filename)
            return False
        
        header = ''
        i_mean, i_variance, i_avgESS, i_psfr = -1,-1,-1,-1
        
        study = min_max + '-' + alig_cons + '-' + std_covar
        self.dic_pstat[study] = {}
        
        for line in f:
            mat = line.rstrip().split('\t')
            
            if mat[0] == 'Parameter':
                if i_psfr == -1:
                    header = deepcopy(mat)
                    
                    for i in range(len(header)):
                        if header[i] == 'Mean':
                            i_mean = i
                        elif header[i] == 'Variance':
                            i_variance = i
                        elif header[i] == 'avgESS':
                            i_avgESS = i
                        elif header[i] == 'PSRF':
                            i_psfr = i
                    
            else:
                for i in range(len(mat)):
                    try:
                        mat[i] = float(mat[i])
                    except:
                        pass
                    
                param = mat[0]
                error = False
                self.dic_pstat[study][param] = {}
                try:
                    self.dic_pstat[study][param]['Mean'] = mat[i_mean]
                except:
                    error = True
                    self.dic_pstat[study][param]['Mean'] = 1
                    
                try:
                    self.dic_pstat[study][param]['Variance'] = mat[i_variance]
                except:
                    error = True
                    self.dic_pstat[study][param]['Variance'] = 1
                    
                try:
                    self.dic_pstat[study][param]['avgESS'] = mat[i_avgESS]
                except:
                    error = True
                    self.dic_pstat[study][param]['avgESS'] = 1
                    
                try:
                    self.dic_pstat[study][param]['PSRF'] = mat[i_psfr]
                except:
                    error = True
                    self.dic_pstat[study][param]['PSRF'] = 1
                    
                if error:
                    print "Error reading mr.bayes files study %s, param %s"%(study,param)
       
        return True

    def read_lstat(self, filename, min_max, alig_cons, std_covar):
        '''
        run    arithmetic_mean    harmonic_mean    values_discarded
        1    -1.477228e+04    -1.494430e+04    no
        2    -1.476670e+04    -1.494118e+04    no
        all    -1.476739e+04    -1.494365e+04    no
        '''

        filename += '.lstat'
        try:
            f = open(filename)
        except:
            print 'could not read %s'%(filename)
            return False
        
        header = ''
        i_arit, i_harmon, i_disc = -1,-1,-1
        
        study = min_max + '-' + alig_cons + '-' + std_covar
        self.dic_ll[study] = {}
        
        ''' run arithmetic_mean    harmonic_mean    values_discarded '''
        for line in f:
            mat = line.rstrip().split('\t')
            
            if mat[0] == 'run':
                if i_disc == -1:
                    header = deepcopy(mat)
                    
                    for i in range(len(header)):
                        if header[i] == 'arithmetic_mean':
                            i_arit = i
                        elif header[i] == 'harmonic_mean':
                            i_harmon = i
                        elif header[i] == 'values_discarded':
                            i_disc = i
                    
            else:
                for i in range(len(mat)):
                    try:
                        if i in [1,2]:
                            mat[i] = float(mat[i])
                    except:
                        pass
                    
                param = mat[0]
                self.dic_ll[study][param] = {}
                self.dic_ll[study][param]['arithmetic_mean'] = mat[i_arit]
                self.dic_ll[study][param]['harmonic_mean'] = mat[i_harmon]
                self.dic_ll[study][param]['values_discarded'] = mat[i_disc]
        
        return True
        
    def stat_files_to_dic(self):

        for min_max in self.min_maxs:
            for alig_cons in self.alig_conss:
                for std_covar in self.std_covars: 
            
                    try:
                        filename = self.mb_filenames[min_max][alig_cons][std_covar]
                    except:
                        continue
                    
                    if os.path.isfile(filename):
                        if not self.read_pstat(filename, min_max, alig_cons, std_covar):
                            return False
                        if not self.read_lstat(filename, min_max, alig_cons, std_covar):
                            return False
                    else:
                        print('>> problems: %s'%(filename))
                        return False
                    
        return True
            
    
    def pstat_critic(self):        
        lista = sorted(self.dic_pstat.keys())
        
        for study in lista:
            lista2 = sorted(self.dic_pstat[study].keys())
            
            print study
            
            for param in lista2:
                try:
                    if self.dic_pstat[study][param]['avgESS'] < 100 or abs(self.dic_pstat[study][param]['PSRF']-1.) > 0.1:
                        print '\t %s = %f (%f) avgESS=%f  PSRF=%f' %(param, self.dic_pstat[study][param]['Mean'], np.sqrt(self.dic_pstat[study][param]['Variance']), \
                                       self.dic_pstat[study][param]['avgESS'], self.dic_pstat[study][param]['PSRF']),
                                       
                        if self.dic_pstat[study][param]['avgESS'] < 100: print ' - problems with avgESS ',             
                        if abs(self.dic_pstat[study][param]['PSRF']-1.) > 0.1: print ' - problems with PSRF ',
                        print ''            
                    else:
                        print '\t %s = %f (%f) avgESS=%f  PSRF=%f' %(param, self.dic_pstat[study][param]['Mean'], np.sqrt(self.dic_pstat[study][param]['Variance']), \
                                       self.dic_pstat[study][param]['avgESS'], self.dic_pstat[study][param]['PSRF'])
                except:
                    print '\t - param not found ', param, self.dic_pstat[study][param]
                
            print '---------------------------'
            
    def likelihood_ratio_test(self, desk, min_max, alig_cons):        
        try:
            H0Stand = min_max + '-' + alig_cons + '-' + self.std_covars[0]
            HaCovar = min_max + '-' + alig_cons + '-' + self.std_covars[1]
            
            if desk.LnL_HMean:
                lnlStd = -hmean(-np.array(self.dic_detaild_params[H0Stand]["LnL"]))
                lnlCov = -hmean(-np.array(self.dic_detaild_params[HaCovar]["LnL"]))
            else:
                lnlStd = np.mean(self.dic_detaild_params[H0Stand]["LnL"])
                lnlCov = np.mean(self.dic_detaild_params[HaCovar]["LnL"])
            
            chi2_calc = -2*(lnlStd - lnlCov) 

            ''' on-off, off -on ''' 
            df = 2  
            chi2_tab = chi2.ppf([0.95], df)
            
            if chi2_calc > chi2_tab:
                stri = "The two models are statistically different with 95% of confidence."
            else:
                stri = "We cannot discard the null hypothesis, we cannot affirm that the two models are different"
                
            return lnlCov, lnlStd, chi2_calc, chi2_tab, stri
            
        except:
            print "Problems with LnL params in LRT."
            return None,None,None,None,None
        

    def write_stat_files(self):  
        for alig_cons in self.alig_conss:
            for min_max in self.min_maxs:
                
                try:
                    _ = self.mb_filenames[min_max][alig_cons]
                except:
                    continue
                
                study =  min_max + '-' + alig_cons
                study_standard = study+'-standard'
                study_covarion = study+'-covarion'
                
                stri = study + '\n'            
                stri += 'Parameter\tConsensus\t\tCovarion\t\tConsensus\t\tCovarion'
                stri += '\n \tmean\tsdv\tmean\tsdv\tavgESS\tPSRF\tavgESS\tPSRF'

                '''---- pstat ----------------------------------------------------'''
                ''' covarion has the on-off parameters '''
                try:
                    params = sorted(self.dic_pstat[study_covarion].keys())
                except:
                    print "dic_pstat does'n has key:", study_covarion
                    continue
                
                for param in params:
                    try:
                        if param == 's(off->on)' or param == 's(on->off)':
                            saux = '\n %s \t  \t '%(param)
                        else:
                            saux = '\n %s \t %e \t %e'%(param, self.dic_pstat[study_standard][param]['Mean'], np.sqrt(self.dic_pstat[study_standard][param]['Variance']) )
                        stri += saux
                    except:
                        pass
                    
                    try:
                        saux = ' \t %e \t %e'%(self.dic_pstat[study_covarion][param]['Mean'], np.sqrt(self.dic_pstat[study_covarion][param]['Variance']) )
                        stri += saux
                    except:
                        pass
                        
                    try:
                        if param == 's(off->on)' or param == 's(on->off)':
                            saux = ' \t  \t '
                            stri += saux
                        else:
                            saux = ' \t %e \t %e'%(self.dic_pstat[study_standard][param]['avgESS'], np.sqrt(self.dic_pstat[study_standard][param]['PSRF']) )
                            stri += saux
                    except:
                        pass
        
                    try:
                        saux = ' \t %e \t %e'%(self.dic_pstat[study_covarion][param]['avgESS'], np.sqrt(self.dic_pstat[study_covarion][param]['PSRF']) )
                        stri += saux
                    except:
                        pass
                    
                '''---- lstat ----------------------------------------------------'''                    
                stri += '\n\n'         

                stri += 'run\tConsensus\t\tCovarion\t\tConsensus\tCovarion'
                stri += '\n \tarithmetic_mean\tharmonic_mean\tarithmetic_mean\tharmonic_mean\tvalues_discarded\tvalues_discarded'

                
                ''' covarion has the on-off parameters '''
                params = sorted(self.dic_ll[study_covarion].keys())
                
                for param in params:
                    saux = '\n %s \t %e \t %e'%(param, self.dic_ll[study_standard][param]['arithmetic_mean'], \
                                                       self.dic_ll[study_standard][param]['harmonic_mean'] )
                    stri += saux
                    
                    saux = '\t%e \t %e'%(self.dic_ll[study_covarion][param]['arithmetic_mean'], \
                                         self.dic_ll[study_covarion][param]['harmonic_mean'] )
                    stri += saux
                        
                        
                    saux = '\t %s \t %s'%(self.dic_ll[study_standard][param]['values_discarded'], \
                                          self.dic_ll[study_covarion][param]['values_discarded'] )
                    stri += saux



                print stri
                
                
                filename = self.desk.mrBayes_path + "mb_params_%s_%s_%s.csv"%(self.desk.organism, self.desk.gene_title, study)
                self.write_file(filename, stri)


    def save_params(self, params):
        stri = "%s \t %i \t %5.3f \t %5.3f \t %5.3f \t %5.3f  \t %5.3f \t %5.3f \t %s \t %5.3f \t %5.3f"
        
        for par in params:
            print "Param: %s"%(par)
            dic = self.summary[par] 
            lista = dic.keys()
            lista = sorted(lista)
            
            print "study".ljust(30) + " \t N   \t mu         \t sigma    \t SE          \t minimum  \t maximum  \t median   \t conf.interval 95% \t avgESS \t PSRF"
            
        
            for min_max in self.min_maxs:
                for alig_cons in self.alig_conss:
                    for std_covar in self.std_covars: 
                        try:       
                            study = min_max + '-' + alig_cons + '-' + std_covar
                            study30 = study.ljust(30)
                            ci = "[%5.3f, %5.3f]"%(dic[study]["q025"], dic[study]["q975"])
                            print stri%(study30, dic[study]["N"], dic[study]["mu"], dic[study]["sigma"], dic[study]["SE"],\
                                        dic[study]["min"], dic[study]["max"], dic[study]["median"], ci,
                                        dic[study]["avgESS"], dic[study]["PSRF"] )
                        except:
                            pass
                
            print "\n-------------------------------\n"
                 
                 
        
        self.pstat_critic()
        self.write_stat_files()
        
                                    
    def write_file(self, filename, stri):
        try:
            f = open(filename, 'w')
            f.write(stri)
            f.flush()
            
            print "File '%s' saved."%(filename)
            ret = True
            
        except:
            print 'Could not write in ' + filename
            ret = False
            
        finally:
            f.close()
            
        return ret

    def show_lnl(self, desk, plt, ax, study, sAvgESS, sPSRF, cntLine, cntCols, min_max, alig_cons, stri_LRT):
        xLabel = "log10(generations)"
        yLabel = "lnl"
        
        mu = self.summary["LnL"][study]["mu"]
        sigma = self.summary["LnL"][study]["sigma"]
        SE = self.summary["LnL"][study]["SE"]
        hmu = self.summary["LnL"][study]["hmu"]
        hSE = self.summary["LnL"][study]["hSE"]
        q2 = self.summary["LnL"][study]["median"]
        
        
        alpha = 0.05
        limInf = round(hmu - norm.ppf(1.-alpha/2.) * hSE,3)
        limSup = round(hmu + norm.ppf(1.-alpha/2.) * hSE,3)
                    
        title = study + "\nmu=%.1f se=%.1f med=%.1f\nhmu=%.1f hse=%.1f CI=[%.1f, %.1f]\nTL %s  %s"%\
                     (mu, SE, q2, hmu, hSE, limInf, limSup, sAvgESS, sPSRF)

        ax.set_title(title, fontsize=9)

        gens = np.log10(self.dic_detaild_params[study]["Gen"])

        vals = self.dic_detaild_params[study]["LnL"]
        y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)
        
        minGen = gens[0]
        maxGen = gens[-1]
        
        plt.plot(gens, vals, 'k-', lw=2, color='blue')
        plt.tick_params(labelsize=9)
        plt.xlim(minGen,maxGen)
        plt.ylim(np.min(vals),np.max(vals))

        #if cntLine == 0:
        plt.xlabel(xLabel, fontsize=9)
        
        if cntLine == 1:
            study2 = min_max + '-' + alig_cons
            lnlCov, lnlStd, chi2Calc, chi2Table, stri = self.likelihood_ratio_test(desk, min_max, alig_cons)
            
                    
            if alig_cons == "aligned":
                stri_LRT += min_max + "\t"
                
            if chi2Calc > chi2Table:
                stri_LRT += "Ha\t"
            else:
                stri_LRT += "H0\t"
                
            if alig_cons == "consensus":
                stri_LRT += "\n"     
                        
                                    
            stri1 = "%s: %s\n%s %6.1f \n%s %6.1f\ncalc chi2=%5.3f chi2(.95,df=2)=%5.3f\n"%\
                  (study2, stri, \
                   "H0=covarion", lnlCov,\
                   "Ha=standard", lnlStd, \
                   chi2Calc, chi2Table)
            print stri1
            
            if chi2Calc > chi2Table: 
                self.test_LRT = "reject H0\nchi2Calc=%5.3f >\nchi2(.95,df=2)=%5.3f"%(chi2Calc, chi2Table[0])
                self.color_LRT = "green"
            else: 
                self.test_LRT = "accept H0\nchi2Calc=%5.3f <=\nchi2(.95,df=2)=%5.3f"%(chi2Calc, chi2Table[0])
                self.color_LRT = "black"
          
            '''
            left = self.dic_lnl[study]["Gen"][12]
            delta = (np.max(self.dic_lnl[study]["LnL"]) - np.min(self.dic_lnl[study]["LnL"])) * .8
            top = np.min(self.dic_lnl[study]["LnL"]) + delta
            ax.text(left, top, self.stri_LRT, color=col2)
            '''
            
        else:
            lnlCov, lnlStd, chi2Calc, chi2Table, stri = -1,-1,-1,-1,""
                                      
        if cntCols == 0:
            plt.ylabel(yLabel, fontsize=9)


        logsup = mu+2*sigma
        loginf = mu-2*sigma

        ax.plot([minGen, maxGen], [hmu, hmu], 'k-', lw=2, color='green')
        ax.plot([minGen, maxGen], [mu, mu], 'k-', lw=2, color='black')
        ax.plot([minGen, maxGen], [logsup, logsup], 'k-', lw=2, color='red')
        ax.plot([minGen, maxGen], [loginf, loginf], 'k-', lw=2, color='red')

        return stri_LRT


    def show_distrib(self, plt, ax, par, study, sAvgESS, sPSRF, ticks, colors, iPar, mini, maxi, cntCols):

        yLabel = "freq(%s)"%par
        mu = self.summary[par][study]["mu"]
        # sigma = self.summary[par][study]["sigma"]
        q025 = self.summary[par][study]["q025"]
        q975 = self.summary[par][study]["q975"]
        q2 = self.summary[par][study]["median"]
        
        vals = self.dic_detaild_params[study][par]

        
        ax.set_title(study + "\nmu=%5.3f med=%5.3f CI=[%5.3f, %5.3f]\n %s  %s"%\
                     (mu, q2, q025, q975, sAvgESS, sPSRF), fontsize=10)

        ax.set_xticks(ticks)
        ax.xaxis.set_ticklabels(ticks)

        # n, bins, patches = ax.hist(vals,  bins=100, color=colors[iPar])  # , facecolor="lightblue"
        n, _, _ = ax.hist(vals,  bins=100, color=colors[iPar])  # , facecolor="lightblue"
        plt.xlim(mini,maxi)

        '''
        if cntLine == 1:
            plt.xlabel(xLabel, fontsize=10)
        '''
        if cntCols == 0:
            plt.ylabel(yLabel, fontsize=10)

                            
        ax.plot([mu, mu], [0, max(n)], 'k-', lw=2, color='black')
        ax.plot([q025, q025], [0, max(n)], 'k-', lw=2, color='red')
        ax.plot([q975, q975], [0, max(n)], 'k-', lw=2, color='red')
        ax.plot([q2, q2], [0, max(n)], 'k-', lw=2, color='yellow')


    def show_conf_interval(self, min_max, alig_cons, par, plt, numLines,numCols, cntCols, ticks, ticks0, mini, maxi, listData):
        ax = plt.subplot2grid((numLines,numCols), (11, cntCols), rowspan=2, colspan=10)
        
        try:
            if cntCols == 0:
                ax.yaxis.set_ticklabels(["","","diff","",""])
            else:
                ax.yaxis.set_ticklabels([])
            plt.ylim(-1,1)
            
            #for x in ticks:
            #    ax.plot([x,x], [-1,1], color="black",linestyle="dotted")
    
            ttest = ttest_ind(listData[0], listData[1], equal_var=False)
            pvalue = ttest[1]
            
            alpha = .05
            if pvalue >= alpha:
                stri = "t-test with identical means, p-value=%.2e"%pvalue
            else:
                stri = "t-test with different means, p-value=%.2e"%pvalue
                
            stri = "Confidence Interval after burnin %.2f\n"%(self.desk.burnin) + stri
            ax.set_title(stri, fontsize=10)
            
            y = 0
       
            study3 = min_max + '-' + alig_cons + '-' + "standard"
            study4 = min_max + '-' + alig_cons + '-' + "covarion"
            
            diff = round(self.summary[par][study4]["mu"] - self.summary[par][study3]["mu"],3)
            N1 = len(self.dic_detaild_params[study3][par])
            var1 = np.var(self.dic_detaild_params[study3][par])
            N2 = len(self.dic_detaild_params[study4][par])
            var2 = np.var(self.dic_detaild_params[study4][par])
            
            SE = np.sqrt(var1/N1 + var2/N2)

            lim1 = round(diff - norm.ppf(1.-alpha/2.) * SE,3)
            lim2 = round(diff + norm.ppf(1.-alpha/2.) * SE,3)
                                    
            ax.plot([lim1,lim2], [y,y], color="red",lw=2)
            ax.plot([lim1,lim1], [y-0.05,y+0.05], color="red",lw=2)
            ax.plot([lim2,lim2], [y-0.05,y+0.05], color="red",lw=2)
            ax.plot([diff,diff], [y-0.05,y+0.05], color="red",lw=2)    
            
            if lim1 > lim2:
                ticks = [lim2,diff,lim1]
            else:
                ticks = [lim1,diff,lim2]
                
            if par == "LnL":
                ticks = [ round(x,1) for x in ticks]

            if ticks[0] > 0: limInf = ticks[0] * .99
            elif ticks[0] < 0: limInf = ticks[0] * 1.01
            
            if ticks[1] < 0: limSup = ticks[2] * .99
            elif ticks[1] > 0: limSup = ticks[2] * 1.01
            
            plt.xlim(limInf,limSup)
            ax.set_xticks(ticks)
            ax.xaxis.set_ticklabels(ticks)
        except:
            print "error in confidence interval"    
       
   
    def show_LRT(self, alig_cons, par, plt, numLines,numCols, cntCols):
        ax = plt.subplot2grid((numLines,numCols), (11, cntCols), rowspan=2, colspan=10)
        
        ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticklabels([])
        plt.ylim(0,1)
                     
        ax.plot([-1,0,1], [0,0,0], color="black")
        
        ax.text(-.8, .3, self.test_LRT , color=self.color_LRT)
        

       
    def calc_ttest_bonferroni(self, data1, data2, pValueLimit):

        ttest = ttest_ind(data1, data2, equal_var=False)
        pvalue = ttest[1]
        
        if pvalue >= pValueLimit:
            return ["Identical", pvalue]
        
        return ["Different", pvalue]
            

            
    def show_qqplot2(self, min_max, alig_cons, par, plt, numLines,numCols, cntCols, listData):
        
        for std_covar3 in ["standard","covarion"]:
            if std_covar3 == "standard":
                x = listData[0]
                ax = plt.subplot2grid((numLines,numCols), (14, cntCols), rowspan=2, colspan=4)
            else:
                x = listData[1]
                ax = plt.subplot2grid((numLines,numCols), (14, cntCols+6), rowspan=2, colspan=4)
            
            study3 = min_max + '-' + alig_cons + '-' + std_covar3
            mu = self.summary[par][study3]["mu"] 
            sigma = self.summary[par][study3]["sigma"]
            # randNorm = np.random.normal(loc=mu, scale=sigma, size=100)   
            ''' Generates a probability plot of sample data against the quantiles of a specified theoretical distribution 
                http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.probplot.html '''
            probplot((x-mu)/sigma, dist="norm", plot=plt)
            ax.set_title("")
            plt.xlabel("")
            plt.ylabel("")
            plt.tick_params(labelsize=10)
            
            if cntCols != 0 or std_covar3 == "covarion":
                ax.yaxis.set_ticklabels([])
            
    def write_newick(self, genetree, treeFilenameNwk):
        try:
           
            dic = {}
            for node in genetree.get_leaves():
                if node.name in dic.keys():
                    dic[node.name] += 1 
                    node.name = str(dic[node.name]) + "." + node.name
                else:
                    dic[node.name] = 1
                    node.name = "1." + node.name
                    
        
            genetree.write(outfile=treeFilenameNwk,format=4)
            ret = True
        except:
            print "could not write newick format"
            ret = False
            
        return ret
                
    def prepare_mr_bayes_data(self):
        desk = self.desk

        if self.desk.aligned_consensus_var.get() == "":
            print "Please confirm: consensus or aligned"
            return False, None
        
        if self.desk.standard_covarion_var.get() == "":
            print "Please confirm: standard or covarion"
            return False, None   
        
        ''' Drosophila_maxmer_Gene_Adh_100L_cutoff7_aligned_covarion.nxs '''
        if desk.standard_covarion_var.get() == "standard":
            if desk.minmax == "mincut":
                if desk.aligned_consensus == "aligned":
                    self.treeFilename = desk.mrBayes_fn01
                else:
                    self.treeFilename = desk.mrBayes_fn02
            else:
                if desk.aligned_consensus == "aligned":
                    self.treeFilename = desk.mrBayes_fn03
                else:
                    self.treeFilename = desk.mrBayes_fn04
        else:
            if desk.minmax == "mincut":
                if desk.aligned_consensus == "aligned":
                    self.treeFilename = desk.mrBayes_fn01_cov
                else:
                    self.treeFilename = desk.mrBayes_fn02_cov
            else:
                if desk.aligned_consensus == "aligned":
                    self.treeFilename = desk.mrBayes_fn03_cov
                else:
                    self.treeFilename = desk.mrBayes_fn04_cov
            
        print self.treeFilename

        maxLnL, maxGen = self.find_max_LnL(self.treeFilename)

        if not maxGen:
            return False, None
            
        print "Found generarion=%i with maxLnl=%e"%(maxGen, maxLnL)

        self.replacements = []
        stri1 = self.desk.mrBayes_replace01_var.get()
        stri2 = self.desk.mrBayes_replace02_var.get()
        stri3 = self.desk.mrBayes_replace03_var.get()
        stri4 = self.desk.mrBayes_replace04_var.get()
        
        for stri in [stri1, stri2, stri3, stri4]:
            if stri != "":
                pos = stri.find(":")
                if pos >= 0:
                    per = stri[pos+1:]
                    stri = stri[:pos]
                    words = stri.split(",")
                    self.replacements.append([words, per])


        for reps in self.replacements:
            for rep in reps[0]:
                print "replace %s per %s"%(rep, reps[1])

        genetree =  self.show_tree(self.treeFilename, maxGen, show=False)
        
        if not genetree:
            return False, None

        self.killed_stri_list = self.desk.string_kill_var.get().split(",")
       
        for i in range(len(self.killed_stri_list)):
            self.killed_stri_list[i] = self.killed_stri_list[i].strip()
            
            
        if self.desk.colaps_tree_var.get():

            nodes = genetree.get_leaves()
            for node in nodes:
                for killed_stri in self.killed_stri_list:
                    node.name = node.name.replace(killed_stri, "")
     
            # L1 = len(nodes)
                     
            while True:
                restart = False
                
                for searchNode in genetree.get_leaves():
                    try:
                        ancestor = searchNode.get_ancestors()[0]
                    except:
                        continue
                    
            
                    leaves = ancestor.get_children()
                    
                    if leaves[0] == searchNode:
                        otherNode = leaves[1]
                    else:
                        otherNode = leaves[0]
                        
                    if not otherNode.is_leaf():
                        continue
    
                        
                    if searchNode.name != otherNode.name:
                        continue
                
    
                    ''' link the first leaf in the pre-ancestor '''
                    # a2.add_child(searchNode)
            
                    name = searchNode.name
                    ''' delete nodes '''
                    searchNode.delete(prevent_nondicotomic=False)
                    otherNode.delete(prevent_nondicotomic=False)
                    ''' ancestor became a node  '''
                    ancestor.name = name
                    
                    # print genetree
                    # print "-----------------"
                    restart = True
            
                if not restart:
                    break
            
            
            # print "---- end ----"

            nodes = genetree.get_leaves()
            # L2 = len(nodes)
            # print "Beginig: num of nodes=%i"%(L1)
            # print "End:     num of nodes=%i"%(L2)
        
        return True, genetree



    def build_see_tree(self):
        ret, genetree = self.prepare_mr_bayes_data()
        
        if not ret:
            return

        treeFilenameNwk = self.desk.mrBayes_path + self.treeFilename.replace(".nxs",".nwk")
        if not self.write_newick(genetree, treeFilenameNwk):
            return
        
        try:
            ts = ete2.TreeStyle()
            ts.show_leaf_name = True
            ts.show_branch_length = True
            ts.show_branch_support = True
            genetree.show(tree_style=ts)
        except:
            if self.treeFilename:
                try:
                    tree = Phylo.read(self.treeFilenameNwk, "newick")
                    Phylo.draw(tree)
                except:
                    print "Could not print in ETE2 neither in Phylo"    
                    print genetree
        
        

    def build_fig_tree(self):
        ret, genetree = self.prepare_mr_bayes_data()
        
        if not ret:
            return
        
        treeFilenameNwk = self.desk.mrBayes_path + self.treeFilename.replace(".nxs",".nwk")
        if not self.write_newick(genetree, treeFilenameNwk):
            return
        
        if not os.path.isfile(treeFilenameNwk):
            print('Could not find: %s'%(treeFilenameNwk))
            return False

        self.desk.call_external_program(self.desk.figtree_path, treeFilenameNwk)

    def save_file(self, fname_path, stri):
        try:
            f = open(fname_path,'w')
            f.write(stri)
            print('Saved %s'%(fname_path))
            ret = True
        except:
            print('Could not write in %s'%(fname_path))
            ret = False
        finally:
            f.close
            
        return ret
            
    def search_correct_tree(self, f, fname_open):
        print 'starting %s'%(fname_open)
        ntree = 1
        first_time = True
        
        count_header = 0
        stri_tree = ''
        stri_back = ''
        list_errors = []
        
        for line in f:
            if count_header < self.num_lines_header_tree:
                stri_tree += line
                stri_back += line
                
                count_header += 1
                continue
    
            pos = line.find(self.mb_clock)
            
            if pos >= 0:
                if first_time: first_time = False; ntree = 1;
                
                line_tree = line[pos + len(self.mb_clock)-1:]    
                

                if ntree%50 == 0:
                    print ntree,
                    if ntree%1200 == 0:
                        print ''
   
                           
                try:
                    ''' see if Tree can be instantiated, otherwise error '''
                    _ = ete2.Tree(line_tree)
                    stri_tree += line
                    stri_back += line
                except:
                    ''' on error: only stri_back '''
                    stri_back += '*** ' + line
                    line_begin = line[:pos]
                    generation = line_begin[self.len_start: ]
                    
                    list_errors.append([generation,line])
                    print '**',
                    
            else:
                if first_time:
                    line = line.lower()
                    line = line.replace(self.del_substring, '')
                    pos1 = line.find(self.del_to_right)
     
                    if pos1 < 0:
                        stri_tree += line
                        stri_back += line
    
                    else:
                        ''' all lines end with , and the last with ; '''
                        line = line[:pos1] + line[len(line)-2:]

                        stri_tree += line 
                        stri_back += line
                else:
                    print '\n\n*** end ***', line
                    stri_tree += line
                    stri_back += line
                
            ntree += 1
    
        if line.lower().find('end;') < 0:
            stri_tree += 'End;\n'
            stri_back += 'End;\n'
    
        return stri_tree, stri_back, list_errors
        
        
    def search_correct_p(self, f, fname_open):
        print 'starting %s'%(fname_open)
        num_p = 1
        
        count_header = 0
        stri_p = ''
        stri_back = ''
        list_errors = []
        
        for line in f:
            if count_header < self.num_lines_header_pi:
                stri_p += line
                stri_back += line
                count_header += 1
                continue
    
            if num_p%50 == 0:
                print num_p,
                if num_p%1200 == 0:
                    print ''         
    
            items = line.split('\t')
            
            if len(items) < 15:
                generation = items[0]
                list_errors.append([generation,line])
                stri_back += '*** ' + line
                print '**',
            else:
                ok = True
                
                for item in items:
                    try:
                        ''' see if Tree can be instantiated, otherwise error '''
                        float(item)
                    except:
                        ok = False
                        break
                        
                if ok:
                    stri_p += line
                    stri_back += line
                else:
                    stri_back += '*** ' + line
                    generation = items[0]
                    
                    list_errors.append([generation,line])
                    print '**',
                          
            num_p += 1
    
        return stri_p, stri_back, list_errors
        
        
    def loop_trees(self):
        for fname in self.fname_nxs:
            for run in ['.run1.t', '.run2.t']:
               
                fname_path = self.desk.mrBayes_path+fname+run
                fname_good = fname_path + '.good'
                fname_back = fname_path + '.backup'
                
                ''' if good does not exists, rename, and use it '''
                if os.path.isfile(fname_path) and os.path.isfile(fname_good):
                    print('Analysis already done: %s'%(fname_good))
                    continue
                
                if not os.path.isfile(fname_good):
                    try:
                        os.rename(fname_path, fname_good)
                        print 'File renamed from \n%s \nto \n%s'%(fname_path, fname_good)
                    except:
                        print 'Error while renaming \n%s \nto \n%s'%(fname_path, fname_good)
                        return False
        
               
                try:
                    f = open(fname_good)
                except:
                    print 'could not read %s'%(fname_good)
                    return False
                        
                stri_tree, stri_back, list_errors = self.search_correct_tree(f, fname_good)
                    
                ''' backup file: with *** and all lines '''
                if not self.save_file(fname_back,stri_back):
                    return False
                ''' tree files eliminate 'bad' trees '''
                if not self.save_file(fname_path,stri_tree):
                    return False
                    
                stri_error = ''
                
                if list_errors:
                    print "Errors:"
                    for lerr in list_errors:
                        print '  > ', lerr[0]
                        stri_error += lerr[0] + '\t' + lerr[1]
                    
                if stri_error != '':
                    filename = fname_path+'.err'
        
                    if not self.save_file(filename,stri_error):
                        exit(-1)
            
                print '----------------------\n'
                
        return True
                
    def verify_tree_num(self, fname, run_num, tree_to_find):
        if run_num == 1:
            run = '.run1.t'
        else:
            run = '.run2.t'
            
         
        fname_path = self.desk.mrBayes_path+fname+run
        
        ''' if good does not exists, rename, and use it '''
        if not os.path.isfile(fname_path):
            print('Did not find %s'%(fname_path))
            exit(-1)
       
        try:
            f = open(fname_path)
        except:
            print 'Could not read %s'%(fname_path)
            exit(-1)
                
        return self.consist_tree(fname_path, f, tree_to_find)



    def consist_tree(self, fname_path, f, tree_to_find, only_species=True):
        print 'starting %s'%(fname_path)
        ntree = 1
        first_time = True
        
        count_header = 0
        taxon = {}
        
        for line in f:
            if count_header < self.num_lines_header_tree:
                count_header += 1
                continue
    
            pos = line.find(self.mb_clock)
            
            if pos >= 0:
                if first_time: first_time = False; ntree = 1;
                
                line_tree = int(line[pos + len(self.mb_clock)-1:] )   
                line_begin = int(line[:pos])
                generation = int(line_begin[self.len_start: ].strip())
                
                if ntree%50 == 0:
                    print ntree,
                    if ntree%1200 == 0:
                        print ''
                                        
                if tree_to_find != generation:
                    ntree += 1
                    continue


   
                try:
                    ''' see if Tree can be instantiated, otherwise error '''
                    genetree = ete2.PhyloTree(line_tree)
                    # t = Tree(line_tree)
                except:
                    ''' on error: only stri_back '''
                    print "problems with %s tree:%s"%( generation, line_tree[:60])
                    exit(-1)
                
                nodes = genetree.get_leaves()
                
                # num_taxa = len(taxon.keys())
                '''  Too many entries in translation table. Maximum number of taxon names to translate is 450 '''
                # numTranslates = len(nodes)
                
                # print 'num_taxa =%i'%(num_taxa)
                # print 'numTranslates =%i'%(numTranslates)
                
                for node in nodes:
                    num_taxon = node.name.strip()
                    
                    if num_taxon not in taxon.keys():
                        print 'Error in taxon %s - %s'%(num_taxon)
                        exit(-1)

                    if only_species:
                        species = taxon[num_taxon]
                        
                        if species.find("bogo") > 0:
                            pass
                        
                        ''' species = xxx_xxxx_species'''
                        pos = species.find("_")
                        if pos > 0:
                            species = species[pos+1:]
                            
                            pos = species.find("_")
                            if pos > 0:
                                species = species[pos+1:]

                        node.name = species
                    else:                        
                        node.name = ' (%s) %s'%(num_taxon, taxon[num_taxon])
                    

                genetree.show()
                #recon_tree.render("phylotree.png", w=750)
    
                    
                return True
            
            else:
                ''' from mrbayes
                    tip nodes have already been indexed 0, 1, 2, ..., numTaxa-1
                '''
                if first_time:
                    line = line.strip().split(' ')
                    num_taxon = line[0]
                    ''' get rid from comma '''
                    desc_taxon = line[1][:len(line[1])-1]
                    
                    taxon[num_taxon] = desc_taxon
                else:
                    print 'Did not find tree %s'%tree_to_find
                    return False
                
            ntree += 1
        return True
    
                
    def loop_ps(self):
        for fname in self.fname_nxs:        
            for run in ['.run1.p', '.run2.p']:
               
                fname_path = self.desk.mrBayes_path+fname+run
                fname_good = fname_path + '.good'
                fname_back = fname_path + '.backup'
                
                ''' if good does not exists, rename, and use it '''
                if os.path.isfile(fname_path) and os.path.isfile(fname_good):
                    print('Analysis already done: %s'%(fname_good))
                    continue
                
                if not os.path.isfile(fname_good):
                    try:
                        os.rename(fname_path, fname_good)
                        print 'File renamed from \n%s \nto \n%s'%(fname_path, fname_good)
                    except:
                        print 'Error while renaming \n%s \nto \n%s'%(fname_path, fname_good)
                        return False
        
               
                try:
                    f = open(fname_good)
                except:
                    print 'could not read %s'%(fname_good)
                    return False
                        
                stri_p, stri_back, list_errors = self.search_correct_p(f, fname_good)
                    
        
                ''' backup file: with *** and all lines '''
                if not self.save_file(fname_back,stri_back):
                    return False
                ''' tree files eliminate 'bad' trees '''
                if not self.save_file(fname_path,stri_p):
                    return False
                    
                stri_error = ''
                if list_errors:
                    print "Errors:"             
                    for lerr in list_errors:
                        print lerr[0]
                        stri_error += lerr[0] + '\t' + lerr[1]
                    
                if stri_error != '':
                    filename = fname_path+'.err'
                    
                    if not self.save_file(filename,stri_error):
                        return False
        
                print '----------------------\n'
        
        return True

    def find_max_LnL(self, filename):
        
        fname_path = self.desk.mrBayes_path + filename + '.run1.p'
     
        if not  os.path.isfile(fname_path):
            print('Did not foind %s'%(fname_path))
            return None, None

        try:
            f = open(fname_path)
        except:
            print 'Could not read %s'%(fname_path)
            return None, None
        
        return self.max_LnL(f)

    def max_LnL(self, f):

        count_header = 0
        maxLnL = -float('inf')
        maxGen = -1
        
        # i = 0
        for line in f:
            
            if count_header < self.num_lines_header_pi:
                # print line
                count_header += 1
                continue
    
            items = line.split('\t')
            # print items
            
            generation = int(items[0])
            lnl = np.double(items[1])
            
            if lnl > maxLnL:
                maxLnL = lnl
                maxGen = generation
                #print ">>>>", maxLnL, maxGen
            '''
            if i > 10:
                exit()
            i += 1
            '''

        return maxLnL, maxGen
        
    def show_tree(self, filename, which_generation, only_species=True, show=True):
        filename = self.desk.mrBayes_path + filename
        
        if not os.path.isfile(filename) :
            print('Tree not found: %s'%(filename))
            return None
        
        print('Tree found: %s'%(filename))
    
        run = '.run1.t'
        filename += run
            
        ''' if good does not exists, rename, and use it '''
        if not os.path.isfile(filename):
            print('Did not find %s'%(filename))
            return None
       
        try:
            f = open(filename)
        except:
            print 'Could not read %s'%(filename)
            return None
    
        print 'starting %s'%(filename)
        first_time = True
        
        count_header = 0
        taxon = {}
        
        num_lines_header_tree = 5
        str_start = '   tree gen.'
        len_start = len(str_start) 
        
       
        for line in f:
            if count_header < num_lines_header_tree:
                count_header += 1
                continue
    
            pos = line.find(self.mb_clock)
            
            if pos >= 0:
                if first_time: first_time = False;
                
                line_tree = line[pos + len(self.mb_clock)-1:]    
                line_begin = line[:pos]
                generation = int(line_begin[len_start: ].strip())
    
                if generation != which_generation:
                    continue
    
                print "found generation %i"%(which_generation)
    
                try:
                    ''' see if Tree can be instantiated, otherwise error '''
                    genetree = ete2.PhyloTree(line_tree)
                    # t = Tree(line_tree)
                except:
                    ''' on error: only stri_back '''
                    print "problems with %s tree:%s"%( generation, line_tree[:60])
                    exit(-1)
                
                nodes = genetree.get_leaves()
                
                # num_taxa = len(taxon.keys())
                '''  Too many entries in translation table. Maximum number of taxon names to translate is 450 '''
                # numTranslates = len(nodes)
                
                # print 'num_taxa =%i'%(num_taxa)
                # print 'numTranslates =%i'%(numTranslates)
                
                for node in nodes:
                    num_taxon = node.name.strip()
                    
                    if num_taxon not in taxon.keys():
                        print 'Error in taxon %s - %s'%(num_taxon)
                        exit(-1)
    
                    if only_species:
                        species = taxon[num_taxon]
                        
                        out = False
                        for reps in self.replacements:
                            for rep in reps[0]:
                                # print "replace %s per %s"%(rep, reps[1])
                                pos = species.find(rep.strip())
                                if pos >= 0:
                                    species = reps[1].strip()
                                    out = True
                                    break
                            if out:
                                break
                            
                        ''' species = xxx_xxxx_species'''
                        pos = species.find("_")
                        if pos > 0:
                            species = species[pos+1:]
                            
                            pos = species.find("_")
                            if pos > 0:
                                pos = species.find("_")
                                species = species[pos+1:]
                        
                        
                        node.name = species
                    else:                        
                        node.name = ' (%s) %s'%(num_taxon, taxon[num_taxon])

                if not show:
                    return genetree

                try:
                    genetree.show()
                except:
                    try:
                        # filename = self.desk.mrBayes_path+"tree.nwk"
                        filename = filename.replace(".nxs.run1.t",".nwk")
                        print filename
                        genetree.write(features=genetree.features, outfile=filename,format=0)
                        tree = Phylo.read(filename, "newick")
                        Phylo.draw(tree)
                    except:
                        print "Could not print in ETE2 neither in Phylo"
                #recon_tree.render("phylotree.png", w=750)
            else:
                ''' from mrbayes
                    tip nodes have already been indexed 0, 1, 2, ..., numTaxa-1
                '''
                if first_time:
                    line = line.strip().split(' ')
                    #print "----------------"
                    #print line
                    num_taxon = line[0]
                    ''' get rid from comma '''
                    try:
                        desc_taxon = line[1][:len(line[1])-1]
                        
                        taxon[num_taxon] = desc_taxon
                    except:
                        pass
                else:
                    print 'Did not find generation=%i'%which_generation
                    return None
                
        return genetree
        
    def verify_generations_trees(self):
        
        for fname in self.fname_nxs:
            for run in ['.run1.t', '.run2.t']:
               
                fname_path = self.desk.mrBayes_path+fname+run
                
                ''' if good does not exists, rename, and use it '''
                if not os.path.isfile(fname_path):
                    print('Could not find tree file: %s'%(fname_path))
                    continue
                
              
                try:
                    f = open(fname_path)
                except:
                    print 'could not read %s'%(fname_path)
                    return False
  
                print 'starting %s'%(fname_path)

                for line in f:
                    pos = line.find(self.mb_clock)
                    
                    if pos >= 0:
                        line_begin = line[:pos]
                        generation = line_begin[self.len_start: ].strip()
                        
                        if generation in self.dic_generation.keys():
                            self.dic_generation[generation] += 1
                        else:
                            self.dic_generation[generation] = 1

                
        return True
    
    
    def verify_generations_pis(self):
      
        for fname in self.fname_nxs:
            for run in ['.run1.p', '.run2.p']:
               
                fname_path = self.desk.mrBayes_path+fname+run
                
                ''' if good does not exists, rename, and use it '''
                if not os.path.isfile(fname_path):
                    print('Could not find p file: %s'%(fname_path))
                    continue
                
              
                try:
                    f = open(fname_path)
                except:
                    print 'could not read %s'%(fname_path)
                    return False
  
                print 'starting %s'%(fname_path)
                count_header = 0

                for line in f:
                    if count_header < self.num_lines_header_pi:
                        count_header += 1
                        continue
            
                    pos = line.find('\t')

                    if pos >= 0:
                        generation = line[:pos].strip()
                        
                        if generation in self.dic_generation.keys():
                            self.dic_generation[generation] += 1
                        else:
                            self.dic_generation[generation] = 1
                       
                    #ntree += 1
                
        return True
    
    
    def saving_common_generations_trees(self):

        for fname in self.fname_nxs:
            for run in ['.run1.t', '.run2.t']:
               
                fname_path = self.desk.mrBayes_path+fname+run
                
                ''' if good does not exists, rename, and use it '''
                if not os.path.isfile(fname_path):
                    print('Could not find tree file: %s'%(fname_path))
                    return False
                
              
                try:
                    f = open(fname_path)
                except:
                    print 'could not read %s'%(fname_path)
                    return False
  
                print 'Searching count == 4 == 2 trees + 2 pis, for %s'%(fname_path)
                stri_tree = ''
                ntree = 1

                for line in f:

                    if ntree%50 == 0:
                        print ntree,
                        if ntree%1200 == 0:
                            print ''
                                            
                    pos = line.find(self.mb_clock)
                    
                    if pos >= 0:
                        line_begin = line[:pos]
                        generation = line_begin[self.len_start: ].strip()
                        
                        if self.dic_generation[generation] == 4:
                            stri_tree += line
                        else:
                            print '** Skip:', line[:60]
                    else:
                        stri_tree += line
                        
                    ntree += 1
                
                ''' tree file saving '''
                print ''
                if not self.save_file(fname_path,stri_tree):
                    return False
        
        return True
    
    def saving_common_generations_pis(self):
      
        for fname in self.fname_nxs:
            for run in ['.run1.p', '.run2.p']:
               
                fname_path = self.desk.mrBayes_path+fname+run
                
                ''' if good does not exists, rename, and use it '''
                if not os.path.isfile(fname_path):
                    print('Could not find p file: %s'%(fname_path))
                    return False
                
              
                try:
                    f = open(fname_path)
                except:
                    print 'could not read %s'%(fname_path)
                    return False
  
                print 'Searching count == 4 2 trees and 2 pis for %s'%(fname_path)
                count_header = 0
                stri_pi = ''
                npi = 1

                for line in f:
                    if count_header < self.num_lines_header_pi:
                        stri_pi += line
                        count_header += 1
                        continue


                    if npi%50 == 0:
                        print npi,
                        if npi%1200 == 0:
                            print ''
                                        
                    pos = line.find('\t')

                    if pos >= 0:
                        generation = line[:pos]
                        
                        if self.dic_generation[generation] == 4:
                            stri_pi += line
                        else:
                            print '** Skip:', line[:60]
                    else:
                        stri_pi += line
                        
                    npi += 1
                
                ''' p file saving '''
                print ''
                if not self.save_file(fname_path,stri_pi):
                    return False
                             
        return True
    
    

    def remaining_branches(self, node, get_rid):
        return node.get_leaves()[0]
        '''
        let_them = []
        
        for node in node.get_leaves():
            if node not in get_rid:
                let_them.append(node)
                
        return let_them    
        '''
    
    
    def ttest_summary(self, par):
        num_comparisons = 12
        pValue = .05
        pValueLimit = pValue / num_comparisons

        stri = "#comparisons=%i, Bonferrroni p-value=%.2e \n"%(num_comparisons, pValueLimit)
        stri += "study       \t| comparison \t| identical\t| p-value\n"
        stri += "----------------------------------------------------------\n"                                        
        striAux2 = "standard x covarion"
        for min_max in self.min_maxs:
            for alig_cons in self.alig_conss:
                striAux1 = min_max + '-' + alig_cons
                study0 = striAux1 + '-' + self.std_covars[0]
                study1 = striAux1 + '-' + self.std_covars[1]
                
                try:
                    ''' return True identical means, False diff means plux p-value '''
                    ret = self.calc_ttest_bonferroni(self.summary[par][study0]["vals"], \
                                                    self.summary[par][study1]["vals"], pValueLimit)
                    
                    stri += "%s \t| %s \t| %s\t| %.3e\n"%(striAux1, striAux2, ret[0], ret[1])
                except:
                    stri += "%s \n"%(striAux1)
         
                                                
        stri += "----------------------------------------------------------\n"                                        
                                                
        striAux2 = "aligned x consensus"
        for min_max in self.min_maxs:
            for std_covar in self.std_covars:
                striAux1 = min_max + '-' + std_covar
                study0 = min_max + '-' + self.alig_conss[0] + '-' + std_covar
                study1 = min_max + '-' + self.alig_conss[1] + '-' + std_covar
                
                try:
                    ''' return True identical means, False diff means plux p-value '''
                    ret = self.calc_ttest_bonferroni(self.summary[par][study0]["vals"], \
                                                    self.summary[par][study1]["vals"], pValueLimit)
                    
                    stri += "%s \t| %s \t| %s\t| %.3e\n"%(striAux1, striAux2, ret[0], ret[1])
                except:
                    stri += "%s \n"%(striAux1)
                                                   
        stri += "----------------------------------------------------------\n"                                        

        striAux2 = "mincut x maxmer"
        for alig_cons in self.alig_conss:
            for std_covar in self.std_covars:
                striAux1 = alig_cons + '-' + std_covar
                study0 = self.min_maxs[0] + '-' + alig_cons + '-' + std_covar
                study1 = self.min_maxs[1] + '-' + alig_cons + '-' + std_covar
                
                try:
                    ''' return True identical means, False diff means plux p-value '''
                    ret = self.calc_ttest_bonferroni(self.summary[par][study0]["vals"], \
                                                    self.summary[par][study1]["vals"], pValueLimit)
                    
                    stri += "%s \t| %s \t| %s\t| %.3e\n"%(striAux1, striAux2, ret[0], ret[1])
                except:
                    stri += "%s \n"%(striAux1)                              
        stri += "----------------------------------------------------------\n"    
        
        
        return stri
        
    