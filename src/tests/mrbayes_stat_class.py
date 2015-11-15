'''
Created on Jul 29, 2015

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica
'''
import platform, os
from copy import deepcopy
import numpy as np #, scipy.stats
from scipy.stats import f_oneway, chi2, norm, ttest_ind, probplot, hmean
from scipy.stats.mstats import mquantiles

class MrBayes_stat_class(object):
    '''
    classdocs
    '''


    def __init__(self):
        self.min_maxs = ['mincut', 'maxmer']
        self.alig_conss = ['aligned', 'consensus']
        self.std_covars = ['standard', 'covarion']
        
        self.name_standard = 'Drosophila_%s_Gene_Adh_100L_cutoff7_%s.nxs'
        self.name_covarion = 'Drosophila_%s_Gene_Adh_100L_cutoff7_%s_covarion.nxs'
        
        organism = "Drosophila"
        seqType = "Gene"
        gene_title = "Adh"
        cutoffLength = 100
        cutoffNumSeq = 7
        self.name_standard = ('%s_%s_%s_%s_%iL_cutoff%i_%s.nxs') %\
             (organism, "%s", seqType, gene_title, cutoffLength, cutoffNumSeq,"%s")
        self.name_covarion = ('%s_%s_%s_%s_%iL_cutoff%i_%s_covarion.nxs') %\
             (organism, "%s", seqType, gene_title, cutoffLength, cutoffNumSeq,"%s")

        
        print self.name_standard
        print self.name_covarion
        #self.name_covarion = 'Drosophila_%s_Gene_Adh_%s_covarion.nxs'
        
        isWindows = (platform.system() == 'Windows') 
        
        if isWindows:
            self.path = "c:\\mrbayes\\"
        else:
            self.path = '/home/flavio/mrbayes/'
        
        self.dic_pstat = {}
        self.dic_ll = {}
        
        '''
        study = [min_max, alig_cons, std_covar]
        
        self.dic_pstat[study]
        
        '''
    def read_runPs(self):
        self.dic_detaild_params = {}
        self.dic_lnl = {}
        
        for alig_cons in self.alig_conss:
            for min_max in self.min_maxs:
                for std_covar in self.std_covars:        
                    if std_covar == 'standard':
                        filename = self.path + self.name_standard%(min_max, alig_cons)
                    else:
                        filename = self.path + self.name_covarion%(min_max, alig_cons)
                    
                    if os.path.isfile(filename):
             
                        #if filename == self.path + "Drosophila_mincut_Gene_Adh_100L_cutoff7_aligned_covarion.nxs":
                        #    filename = self.path + "Drosophila_mincut_Gene_Adh_100L_cutoff7_aligned_covarion_temp005.nxs"
                                        
                        print('Ok: %s'%(filename))
                        if not self.read_run1p(filename, min_max, alig_cons, std_covar):
                            return False
                            
        return True
                             
    def stat_ppf(self, x):
        x1, x2, x3, x4, x5 = mquantiles(x, prob=[0.025, 0.25, 0.5, 0.75, .975])
        return np.mean(x), np.sqrt(np.var(x)), x1, x2, x3, x4, x5       
           

    def my_ppf(self, x):
        x = sorted(x)
        L = len(x)

        mu = np.mean(x)
        sigma = np.sqrt(np.var(x))

        qls = []
        for alpha in [.025, .25, .5, .75, .975]:
            La = alpha*(L+1)
            diff = La - np.floor(La)
            
            if int(np.floor(La)) < 1:
                qls.append( x[0])
            else:
                floor = int(np.floor(La)) - 1
                
                if floor >= L-1:
                    qls.append(x[L-1])
                else:
                    ceil = int(np.ceil(La)) - 1
                    qls.append(x[floor] + (x[ceil] - x[floor]) * diff)

        return mu, sigma, qls[0], qls[1], qls[2], qls[3], qls[4] 
    
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
         
                                
    def read_run1p(self, filename, min_max, alig_cons, std_covar, burnin=.25):
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

        ini = int((len(lines)-2)*burnin)+1
        
        for par in header:
            self.dic_detaild_params[study][par] = []
        
        for line in lines[ini:]:
            mat = line.rstrip().split('\t')
            if len(mat) != len(header):
                pass
                continue
            
            for i in range(len(mat)):
                par = header[i]
                try:
                    if i==0:
                        self.dic_detaild_params[study][par].append(int(mat[i]))
                    else:
                        self.dic_detaild_params[study][par].append(float(mat[i]))
                except:
                    if i==0:
                        print "found some error converting Generation like integer",str(mat[i])
                    else:
                        print "found some error converting values like float",str(mat[i])
                

        self.dic_lnl[study] = {}
        self.dic_lnl[study]["Gen"] = []
        self.dic_lnl[study]["LnL"] = []
        
        for line in lines[10:]:
            mat = line.rstrip().split('\t')
            if len(mat) != len(header):
                pass
                continue

            try:
                self.dic_lnl[study]["Gen"].append(np.log10(int(mat[0])))
            except:
                print "found some error converting Generation like integer",str(mat[0])
            try:
                self.dic_lnl[study]["LnL"].append(np.log10(-float(mat[1])))
            except:
                print "found some error converting LnL like float",str(mat[1])

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
                self.dic_pstat[study][param] = {}
                self.dic_pstat[study][param]['Mean'] = mat[i_mean]
                self.dic_pstat[study][param]['Variance'] = mat[i_variance]
                self.dic_pstat[study][param]['avgESS'] = mat[i_avgESS]
                self.dic_pstat[study][param]['PSRF'] = mat[i_psfr]

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

        for alig_cons in self.alig_conss:
            for min_max in self.min_maxs:
                for std_covar in self.std_covars:        
                    if std_covar == 'standard':
                        filename = self.path + self.name_standard%(min_max, alig_cons)
                    else:
                        filename = self.path + self.name_covarion%(min_max, alig_cons)
                    
                    if os.path.isfile(filename):
                        # print('Ok: %s'%(filename))
                        #if filename == self.path + "Drosophila_mincut_Gene_Adh_100L_cutoff7_aligned_covarion.nxs":
                        #    filename = self.path + "Drosophila_mincut_Gene_Adh_100L_cutoff7_aligned_covarion_temp005.nxs"
                            
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
            
    def likelihood_ratio_test(self, min_max, alig_cons):        
        try:
            H0Stand = min_max + '-' + alig_cons + '-' + self.std_covars[0]
            lnlStd = np.mean(self.dic_detaild_params[H0Stand]["LnL"])
            
            HaCovar = min_max + '-' + alig_cons + '-' + self.std_covars[1]
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
                study =  min_max + '-' + alig_cons
                study_standard = study+'-standard'
                study_covarion = study+'-covarion'
                
                stri = study + '\n'            
                stri += 'Parameter\tConsensus\t\tCovarion\t\tConsensus\t\tCovarion'
                stri += '\n \tmean\tsdv\tmean\tsdv\tavgESS\tPSRF\tavgESS\tPSRF'

                '''---- pstat ----------------------------------------------------'''
                ''' covarion has the on-off parameters '''
                params = sorted(self.dic_pstat[study_covarion].keys())
                
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
                
                filename = self.path + study + '.csv'
                self.write_file(filename, stri)
                            
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
    
    def show_lnl(self, plt, ax, dic, study, mu, sigma, q2, q025, q975, sAvgESS, sPSRF, vals, cntLine, cntCols, min_max, alig_cons, stri_LRT):
        xLabel = "log10(generations)"
        yLabel = "log10(-lnl)"
        
        ax.set_title(study + "\nmu=%6.0f med=%6.0f CI=[%6.0f, %6.0f]\n TL %s  %s"%\
                     (mu, q2, q025, q975, sAvgESS, sPSRF), \
                     fontsize=10)
        
        vals = [np.log10(-x) for x in vals]

        plt.plot(self.dic_lnl[study]["Gen"], self.dic_lnl[study]["LnL"], 'k-', lw=2, color='blue')

        if cntLine == 0:
            plt.xlabel(xLabel, fontsize=10)
        elif cntLine == 1:
            study2 = min_max + '-' + alig_cons
            lnlCov, lnlStd, chi2Calc, chi2Table, stri = self.likelihood_ratio_test(min_max, alig_cons)
            
                    
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
                stri2 = "reject H0\nchi2Calc=%5.3f >\nchi2(.95,df=2)=%5.3f"%(chi2Calc, chi2Table[0])
                col2 = "green"
            else: 
                stri2 = "accept H0\nchi2Calc=%5.3f <=\nchi2(.95,df=2)=%5.3f"%(chi2Calc, chi2Table[0])
                col2 = "black"
          
            left = self.dic_lnl[study]["Gen"][10]
            top = self.dic_lnl[study]["LnL"][4]
            ax.text(left, top, stri2, color=col2) 
            
        else:
            lnlCov, lnlStd, chi2Calc, chi2Table, stri = -1,-1,-1,-1,""
                                      
        if cntCols == 0:
            plt.ylabel(yLabel, fontsize=10)
        
        minGen = self.dic_lnl[study]["Gen"][0]
        maxGen = self.dic_lnl[study]["Gen"][len(dic["Gen"])-1]
        logmu = np.log10(-mu)
        logsup = np.log10(-mu-sigma)
        loginf = np.log10(-mu+sigma)

        ax.plot([minGen, maxGen], [logmu, logmu], 'k-', lw=2, color='black')
        ax.plot([minGen, maxGen], [logsup, logsup], 'k-', lw=2, color='black',linestyle="dotted")
        ax.plot([minGen, maxGen], [loginf, loginf], 'k-', lw=2, color='black',linestyle="dotted")
        
        return stri_LRT

    
    def show_distrib(self, plt, ax, par, study, mu, sigma, q2, q025, q975, sAvgESS, sPSRF, ticks, vals, colors, iPar, mini, maxi, cntCols):
        # xLabel = par
        yLabel = "freq(%s)"%par
        
        ax.set_title(study + "\nmu=%5.3f med=%5.3f CI=[%5.3f, %5.3f]\n %s  %s"%\
                     (mu, q2, q025, q975, sAvgESS, sPSRF), \
                     fontsize=10)

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
        ax.plot([mu+sigma, mu+sigma], [0, max(n)], 'k-', lw=2, color='black',linestyle="dotted")
        
        if mu-sigma >= 0:
            ax.plot([mu-sigma, mu-sigma], [0, max(n)], 'k-', lw=2, color='black',linestyle="dotted")
        ax.plot([mu+2*sigma, mu+2*sigma], [0, max(n)], 'k-', lw=2, color='red',linestyle="dotted")
        
        if mu-2*sigma >= 0:
            ax.plot([mu-2*sigma, mu-2*sigma], [0, max(n)], 'k-', lw=2, color='red',linestyle="dotted")
            
        ax.plot([q025, q025], [0, max(n)], 'k-', lw=2, color='red')
        ax.plot([q975, q975], [0, max(n)], 'k-', lw=2, color='red')
        ax.plot([q2, q2], [0, max(n)], 'k-', lw=2, color='yellow')


    def show_conf_inteval(self, min_max, alig_cons, par, summary, plt, numLines,numCols, cntCols, ticks, ticks0, mini, maxi, listData):
        ax = plt.subplot2grid((numLines,numCols), (11, cntCols), rowspan=2, colspan=2)
        
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
                stri = "t-test with identical means, p-value=%1.3e"%pvalue
            else:
                stri = "t-test with different means, p-value=%1.3e"%pvalue
                
            stri = "Confidence Interval after burnin 25%\n" + stri
            ax.set_title(stri, fontsize=10)
            
            y = 0
       
            study3 = min_max + '-' + alig_cons + '-' + "standard"
            study4 = min_max + '-' + alig_cons + '-' + "covarion"
            
            diff = round(summary[par][study4]["mu"] - summary[par][study3]["mu"],3)
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
       
   
   
    def show_conf_inteval_3(self, min_max, alig_cons, par, summary, plt, numLines,numCols, cntCols, ticks, ticks0, mini, maxi, listData):
        ax = plt.subplot2grid((numLines,numCols), (11, cntCols), rowspan=2, colspan=2)
        
        try:
            #ax.plot(ticks, ticks0, color="black")
            #plt.xlim(mini,maxi)
            #ax.set_xticks(ticks)
            #ax.xaxis.set_ticklabels(ticks)
            
            if cntCols == 0:
                ax.yaxis.set_ticklabels(["","covar","stand","diff",""])
            else:
                ax.yaxis.set_ticklabels([])
            plt.ylim(-1,1)
            
            #for x in ticks:
            #    ax.plot([x,x], [-1,1], color="black",linestyle="dotted")
    
            ttest = ttest_ind(listData[0], listData[1], equal_var=False)
            pvalue = ttest[1]
            if pvalue >= .05:
                stri = "t-test with identical means, p-value=%1.3e"%pvalue
            else:
                stri = "t-test with different means, p-value=%1.3e"%pvalue
                
            stri = "Confidence Interval after burnin 25%\n" + stri
            ax.set_title(stri, fontsize=10)
            
            for y, stand_covar in [ [-0.5, "covarion"],[0, "standard"],[+0.5, "difference"]]:
       
                if stand_covar != "difference":
                    study3 = min_max + '-' + alig_cons + '-' + stand_covar
                    
                    mu = summary[par][study3]["mu"] 
                    q025 = summary[par][study3]["q025"] 
                    q975 = summary[par][study3]["q975"] 
                
                    ax.plot([q025,q975], [y,y], color="red",lw=2)
                    ax.plot([q025,q025], [y-0.05,y+0.05], color="red",lw=2)
                    ax.plot([q975,q975], [y-0.05,y+0.05], color="red",lw=2)
                    ax.plot([mu,mu], [y-0.05,y+0.05], color="red",lw=2)
                
                else:
                    study3 = min_max + '-' + alig_cons + '-' + "standard"
                    study4 = min_max + '-' + alig_cons + '-' + "covarion"
                    
                    diff = summary[par][study4]["mu"] - summary[par][study3]["mu"]
                    N1 = len(self.dic_detaild_params[study3][par])
                    var1 = np.var(self.dic_detaild_params[study3][par])
                    N2 = len(self.dic_detaild_params[study4][par])
                    var2 = np.var(self.dic_detaild_params[study4][par])
                    
                    SE = np.sqrt(var1/N1 + var2/N2)

                    lim1 = diff - 2 * SE
                    lim2 = diff + 2 * SE
                                            
                    ax.plot([lim1,lim2], [y,y], color="red",lw=2)
                    ax.plot([lim1,lim1], [y-0.05,y+0.05], color="red",lw=2)
                    ax.plot([lim2,lim2], [y-0.05,y+0.05], color="red",lw=2)
                    ax.plot([diff,diff], [y-0.05,y+0.05], color="red",lw=2)    
                      
            
        except:
            print "error in confidence interval"    
       
            
  
            
    def show_qqplot2(self, min_max, alig_cons, par, summary, plt, numLines,numCols, cntCols, listData):
        
        for std_covar3 in ["standard","covarion"]:
            if std_covar3 == "standard":
                x = listData[0]
                ax = plt.subplot2grid((numLines,numCols), (14, cntCols), rowspan=2)
            else:
                x = listData[1]
                ax = plt.subplot2grid((numLines,numCols), (14, cntCols+1), rowspan=2)
            
            study3 = min_max + '-' + alig_cons + '-' + std_covar3
            mu = summary[par][study3]["mu"] 
            sigma = summary[par][study3]["sigma"]
            # randNorm = np.random.normal(loc=mu, scale=sigma, size=100)   
            probplot((x-mu)/sigma, dist="norm", plot=plt)
            ax.set_title("")
            plt.xlabel("")
            plt.ylabel("")
            
            if cntCols != 0 or std_covar3 == "covarion":
                ax.yaxis.set_ticklabels([])
            