#-*- coding: utf-8 -*-
'''
Created on 18/04/2013
Updated on 23/06/2014 - var(JSD) revised
Updated on 23/12/2014 - distance =  sqrt(Jensen�Shannon divergence)
Updated on 22/09/2015 - calc_JS_Distance_simple - for hypothesis test

@author: Flavio Lichtenstein
@local: Unifesp DIS Bioinformatica

@site: cross validate
@url: http://stats.stackexchange.com/questions/29578/jensen-shannon-divergence-calculation-for-3-prob-distributions-is-this-ok

'''
import numpy as np
import BarGraphic as graphPac
import Timing as timePack
import matplotlib.pyplot as plt
import warnings
#from matplotlib.tight_bbox import adjust_bbox

class Class_JSD:
    def __init__(self, desk):
        # self.base = base
        self.desk = desk
        
        self.clearEntropySeqs()
        self.time = timePack.Timing()
        
        self.myPlot = graphPac.Plot()

    def clearEntropySeqs(self):
        self.entropyShannon = [] 
        
    def entropy(self, prob_dist):
        return -sum([p * np.log(p) for p in prob_dist if (p > 0 and p < 1)])
    
    def entropy_multi_line(self, bi_dim_prob_dist):
        h = []
        for prob_dist in bi_dim_prob_dist:
            h.append(-sum([p * np.log(p,self.base) for p in prob_dist if (p > 0 and p < 1)]))
        return h
    

    def calc_JSD(self, listSpecies, list_distMI, list_distSE):
        numSpecies = len(list_distMI)
        roundVal=4
        
        desk = self.desk

        stri = 'species(%s)'%(desk.unit)
        sError = 'species'
        
        dic = {}
        for i in range(numSpecies):
            stri += '\t' + listSpecies[i]
            sError += '\t' + listSpecies[i]
            
        for i in range(numSpecies):
            stri += '\n' + listSpecies[i]
            sError += '\n' + listSpecies[i]

            for j in range(numSpecies):
               
                if i == j:
                    d, dError = 0.0, 0.0
                else:
                    spec1 = listSpecies[i] 
                    spec2 = listSpecies[j] 
                    
                    if i > j:
                        d,dError = dic[spec2 +' x '+ spec1]
                        dic[spec2 +' x '+ spec1] = [d,dError]
                    else:
                        d, dSup, dInf = self.calc_JS_Distance(list_distMI[i],list_distSE[i], list_distMI[j],list_distSE[j])
                        
                        ''' The square root of the Jensen-Shannon divergence is a  
                            metric often referred to as Jensen-Shannon distance '''
                        d    = round(d   ,roundVal)
                        dSup = round(dSup,roundVal)
                        dInf = round(dInf,roundVal)

                        dError = round( abs(dSup - dInf) / 2.,roundVal)
                        dic[spec1+' x '+spec2] = [d,dError]
                
                stri += '\t' + str(d)
                sError += '\t' + str(dError)

        # print '\n',stri,'\n'
        return dic, stri, sError


    
    def calc_JS_Distance(self, distMI1, distSE1, distMI2, distSE2):
        ''' PIj are all equal 
            JSD = H(sum(PIi * pi)) - Sum(PIi H(Pi)) 
        '''
        desk = self.desk
        numOfSEs=desk.numOfSEs
        
        # _ = self.time.start()
        #JSD_mean1 = self.calc_jsd_2dist_wo_numpy(distMI1, distMI2)
        #_ = self.time.finish()
        #t1 = self.time.milli()
        
        
        #_ = self.time.start()
        JS_Dist_mean = self.calc_JS_Dist_2dist(distMI1, distMI2)
        #_ = self.time.finish()
        #t2 = self.time.milli()

        #stri = 'JS_Dist_mean calc_jsd_2dist_wo_numpy %f, %i ms'%(JSD_mean1, t1)
        #self.desk.showmsg_obs(stri)
        #stri = 'JS_Dist_mean jsd %f, %i ms'%(JS_Dist_mean, t2)
        #self.desk.showmsg_obs(stri)
        self.desk.showmsg_obs('.', same_line=True)


        listJS_Dist = []
        #_ = self.time.start()
        listJS_Dist.append(self.calc_JS_Dist_2dist(distMI1+numOfSEs*distSE1,  distMI2+numOfSEs*distSE2))
        listJS_Dist.append(self.calc_JS_Dist_2dist(distMI1+numOfSEs*distSE1,  distMI2-numOfSEs*distSE2))
        listJS_Dist.append(self.calc_JS_Dist_2dist(distMI1-numOfSEs*distSE1,  distMI2+numOfSEs*distSE2))
        listJS_Dist.append(self.calc_JS_Dist_2dist(distMI1-numOfSEs*distSE1,  distMI2-numOfSEs*distSE2))
        #_ = self.time.finish()
        #t4 = self.time.milli()
        
        JS_Dist_max  = max(listJS_Dist)
        JS_Dist_min  = min(listJS_Dist)
        
        #stri = 'JSD max=%f min=%f, %i ms'%(JS_Dist_max, JS_Dist_min, round(t4/4.))
        #self.desk.showmsg_obs(stri)
        
        return JS_Dist_mean, JS_Dist_max, JS_Dist_min

    def calc_JS_Distance_simple(self, distMI1, distMI2):
        JS_Dist_mean = self.calc_JS_Dist_2dist(distMI1, distMI2)
      
        return JS_Dist_mean
    
    
    def calc_JS_Dist_2dist(self, P,Q): #Jensen-shannon divergence
        
        warnings.filterwarnings("ignore", category = RuntimeWarning)
        
        P = np.array(P)
        Q = np.array(Q)
        d1 = P*np.log(2*P/(P+Q))
        d2 = Q*np.log(2*Q/(P+Q))
        d1[np.isnan(d1)] = 0
        d2[np.isnan(d2)] = 0

        return np.sqrt(0.5*np.sum(d1+d2))
        
    ''' must receive normalized data  '''
    def calc_jsd_2dist_wo_numpy(self, P,Q):
        ''' all same weight - here weight = .5 '''
        weight = .5
        l = len(P)
        js_left = np.zeros(l)
        js_right = 0.

        for pd in [P,Q]:
            for k in range(l):
                js_left[k] += pd[k]

            js_right += self.entropy(pd)

        '''    H(Sum (PIj * pj) - Sum(PIj.H(pj) '''
        return self.entropy(js_left*weight)-(js_right*weight)
            

    def calcJSD_sde(self, autocorrSpec0, autocorrSpec1):
        ''' HMI np.mean and sdv for each k from n sequences '''
        HXY = self.entropy_multi_line(autocorrSpec0 + autocorrSpec1)
        HX  = self.entropy_multi_line(autocorrSpec0)
        HY  = self.entropy_multi_line(autocorrSpec1)
                
        sde = np.sqrt(np.std(HXY)**2/len(HXY) + np.std(HX)**2/len(HX) + np.std(HY)**2/len(HY))

        return sde
        
    def save_file(self, filename, stri):
        try:
            f = open(filename, 'w')
            s = str(stri)
            f.write(s)
            f.flush()
            
            ret = True
        except:
            ret = False
        finally:
            f.close()
            
        return ret
    
    '''                    desk, dicMI, self.factor, title, xlabel, ylabel, roundVal=self.roundVal, filename=filename_dist_ma'''
    def plot_JSD_MI(self, desk, dicDist, xlabel, ylabel = 'JSD', onlyRef=False, stay=False):
        factor = desk.factor
        roundVal = desk.roundVal
                
        title = desk.title_jsd + ', bar error=%i*SE'%(desk.numOfSEs)
        i = 0
        lista = sorted(dicDist.keys())
        
        seqX = []
        seqY = []
        seqError = []
        seqSeq = []

        has_negative = False
        for key in lista:
            if onlyRef:
                if key.find("**") < 0: continue
                
            if desk.isLog:
                try:
                    d = np.array(dicDist[key][0])
                    se = np.array(dicDist[key][1])
                    
                    dist = np.round(np.log10(d),roundVal)
               
                    if not has_negative and dist < 0:
                        has_negative = True
                except:
                    dist = 0


                try:
                    distError = np.abs( np.round(np.log10(d)     ,roundVal) - \
                                        np.round(np.log10(d+se)) ,roundVal)
                    
                except:
                    distError = 0
            else:
                d = np.array(dicDist[key][0])
                se = np.array(dicDist[key][1])
                
                if factor == 1:
                    dist      = np.round(d,roundVal)
                    distError = np.round(se,roundVal)
                else:
                    dist      = np.round(d*factor,roundVal)
                    distError = np.round(se*factor,roundVal)
                
            if desk.numOfSEs != 1:
                distError *= desk.numOfSEs
            
            seqX.append(i)
            
            if len(key) <= 22:
                key2 = key
            else:
                elems = key.split(' x ')
                key2 = elems[0][:13].rjust(13, ' ') + '\n' + elems[1][:13].rjust(13, ' ')
            

            seqSeq.append(key2)
            seqY.append(dist)
            seqError.append(distError)
            i += 1

        
        '''
        _, ax = plt.subplots()
        ax.bar(seqX, seqY, yerr=seqError,color='y')
        ax.set_xticklabels( seqSeq )
        ax.set_xticks(seqX)
        '''
        try:
            plt.close("all")
            plt.clf() 
        except:
            pass           
        fig = plt.figure(1, dpi=desk.dpi)  # figsize=(12.0, 8.0), 
        ax = fig.add_subplot(111)
        
        ax.bar(seqX, seqY, yerr=seqError,color='y')
        ax.set_xticklabels(seqSeq)
        ax.set_xticks(seqX)


        mean = np.mean(seqY)
        sd = np.sqrt(np.var(seqY))
        stdY2 = 2 * sd
                        
        seqSup = []
        seqMean = []
        seqInf = []
        for _ in range(len(seqY)):
            seqSup.append(mean+stdY2)
            seqMean.append(mean)
            seqInf.append(mean-stdY2)
        
        plt.plot(seqX, seqSup, color='red')
        plt.plot(seqX, seqMean, color='black')
        plt.plot(seqX, seqInf, color='red')
        
        title += "\n mean(JSD[HMI]) = %.1f and sd(JSD[HMI]) = %.1f"%(mean, sd)
  

        plt.xlim(0, len(seqX) )

        if has_negative: plt.ylim(min(seqY), max(seqY))
       
        fontsize = 22
        plt.title(title, fontsize=fontsize)
        plt.xlabel(xlabel, fontsize=fontsize)
        plt.ylabel(ylabel, fontsize=fontsize)
        
         
        if len(seqX) <= 30: labelsize = 14
        elif len(seqX) <= 60: labelsize = 12
        else: labelsize = 10
        
        plt.xticks(rotation=90)
        plt.subplots_adjust(top=.87, bottom=.22, left=.08, right=.96)        
        # plt.tick_params(axis='both', which='major', labelsize=labelsize)
        plt.tick_params(labelsize=labelsize)

        '''
        if desk.result:
            seqX.append(max(seqX)+1)
            
            ax.plot(seqX, [desk.result["mean"] for _ in range(len(seqX))], color='black')
            ax.plot(seqX, [desk.result["median"] for _ in range(len(seqX))], color='orange')
            ax.plot(seqX, [desk.result["pvalue_inf"] for _ in range(len(seqX))], color='red')
            ax.text(7, desk.result["pvalue_inf"]+5, desk.result["result"], ha='center', fontsize=18, color="blue")
        '''
        pictureName = desk.prefix + desk.sufix + desk.filename_correction
        self.myPlot.print_graph(desk, fig, pictureName=pictureName, frame=desk.tk_root, stay=stay)
  
    def calcVerticalMutualInfo_Dirichlet(self, freqDistrib, numOfLetters=1, showmessage=False):
        ''' -------------------------------------------------------
        --- Calculating entropy
        ------------------------------------------------------- '''
        seqMIijs = []
        seqHis = []
        seqHijs = []
        seqPijs = []
        
        l = len(freqDistrib['A'])

        
        seqPis, seqVars = self.calcShannonVerticalEntropy_Dirichlet(freqDistrib, numOfLetters, showmessage=False)
        
        # meanPis = list(map(np.mean, zip(*seqPis)))
        # stdvPis = list(map(np.std,  zip(*seqPis)))
        
        print '- Calc join Entropies --- len entropies:', len(self.entropyShannon)
        # maxi = 0
        
        for i in range(len(self.entropyShannon) - (numOfLetters-1) ):
            seqHis.append(self.entropyShannon[i])
            if i%50 == 0:
                print i, 
            for j in range(i+1, len(self.entropyShannon)):
                
                h1 = self.entropyShannon[i]
                se = 0
                for var in seqVars[i]:
                    se  += var

                h2 = self.entropyShannon[j]
                for var in seqVars[j]:
                    se  += var


                h12, pijs, var_ij = self.calcShannonDoubleVerticalEntropy_Dirichlet(freqDistrib, i, j, numOfLetters)
                
                se  += var_ij
                
                se = np.sqrt(se /l)
                # print i,j, h1, h2, h12, 'se', se

                seqPijs.append(pijs)
                mi = h1 + h2 - h12

                if mi <> 0:
                    if mi < 1.e-8 or mi - se <= 0:
                        mi = 0.
                        se = 0.
                    '''
                    else:
                        if mi > maxi:
                            maxi = mi
                            print '>>> ',
                            
                        print i, j, 'MI:', mi, 'se', se, '>> h1', h1, 'h2', h2, 'h12', h12
                    '''
                else:
                    se = 0.
                    
                ''' codigo novo como doGraph3D '''
                seqMIijs.append([i, j, mi, se])  
                seqHijs.append([i, j, h1, h2, h12])
                        
                '''
                if h12 <> 0:
                    print count, i, j, mi, '=', h1, '+', h2, '-', h12
                '''
        
        print '\n Pijs len', len(pijs), ' len seqMIijs', len(seqMIijs), 'seqHijs', len(seqHijs)
        # meanPijs = list(map(np.mean, zip(*seqPijs)))
        # stdvPijs = list(map(np.std,  zip(*seqPijs)))


        return seqMIijs, seqHijs, seqVars, seqHis, seqPis, seqPijs


    def calcShannonVerticalEntropy_Dirichlet(self, freqDistrib, numOfLetters, showmessage=False):
        self.clearEntropySeqs()
        
        self.numSequences = len(freqDistrib['A'])
            
        seqPis = []
        seqVars = []       
        for i in range(len(self.numSequences)): 
            ''' mats has frequencies from A, G, T, C '''
            hShan = 0
            # print '...', i, mats

            seq_aux = []
            seq_auxVar = []

            for nuc in ['A','G','T','C']:
                # print mat
                pi = freqDistrib[nuc][i]
                seq_aux.append(pi)

                ''' var(f(x)) = f'(<x>)^2 * Var(X)
                    f(X) = h(p)
                    f'(X) = h'(p) = -(ln pi + 1/ln2)
                '''
    
                if (pi != 0.) and (pi != 1.):
                    seq_auxVar.append(pi * (1-pi) * (np.log(pi,2) + self.basic.inv_ln2) ** 2)
                    hShan -= pi * np.log10(pi)
                else:
                    seq_auxVar.append(0.)
                    
                    # print '<<<', i, mat, pi
                    # print i,mat #',', 

            seqPis.append(seq_aux)
            seqVars.append(seq_auxVar)
            # print 'increment entropy ', m, '  loop ', i
            
            if hShan == 0:
                self.entropyShannon.append(hShan)
            else:
                self.entropyShannon.append(hShan / self.basic.log2)

        # print '--- end counting ---- init hShannon---------'

        if showmessage:
            stri = 'Calculating Shannon'
                
            print '\n------------------------------------------------'         
            print stri + ' Entropy '         
            print '--------------------------------------------------\n'         
            
            for i in range(len(self.entropyShannon)):
                if self.entropyShannon[i] > 0:
                    print i, ' hShan=', self.entropyShannon[i]
                    
        return seqPis, seqVars


    def calcShannonDoubleVerticalEntropy_Dirichlet(self, freqDistrib, posI, posJ, numOfLetters):

        dic = self.countTwoSiteWords(freqDistrib, posI, posJ, numOfLetters) 
        '''
        if posI == 9 and posJ == 371:
            for key in dic.keys():
                print key, dic[key]
                
            print 'yes...'
        '''
        
        seqPijs = []
        hShan = 0 
        var = 0

        for key in dic.keys():
            pi = dic[key][2]
            seqPijs.append(pi)

            if (pi != 0.) and (pi != 1.):
                hShan -= pi * np.log(pi,2)
                
                var += pi * (1-pi) * (np.log(pi) + self.basic.inv_ln2) ** 2

        return hShan, seqPijs, var
    