'''
Created on Jul 29, 2015

@author: flavio
'''
from tests import mrbayes_stat_class as mb
from numpy import sqrt
import matplotlib.pyplot as plt
import numpy as np

mbs = mb.MrBayes_stat_class()

x = list(range(1,11))
print x
mu, sigma, x1, x2, x3, x4, x5 = mbs.stat_ppf(x)
print x
print mu, sigma, x1, x2, x3, x4, x5
print "---"
mu, sigma, x1, x2, x3, x4, x5 = mbs.my_ppf(x)
print x
print mu, sigma, x1, x2, x3, x4, x5
print "---"
mu, sigma, x1, x2, x3, x4, x5 = mbs.my_harmonic_ppf(x)
print x
print mu, sigma, x1, x2, x3, x4, x5
print "---"

exit()

    

if not mbs.read_runPs():
    exit()

summary = {}

params = ["pi(A)","pi(C)","pi(G)","pi(T)"]
params = ["pi(A)","pi(C)"]
params = ["r(A<->C)"]
params = ["pi(A)","pi(C)","s(off->on)","s(on->off)","pinvar"]
params = ["s(off->on)","s(on->off)","pinvar"]
params = ["pi(A)","pi(C)","pi(G)","pi(T)","r(A<->C)","r(A<->G)","r(A<->T)","r(C<->G)","r(C<->T)","r(G<->T)","LnL","LnPr","TL","alpha","s(off->on)","s(on->off)","pinvar"]
params = ["LnL"]
params = ["pi(A)"]
params = ["pi(A)","r(A<->C)","pinvar","LnL"]

''' http://stackoverflow.com/questions/12439588/how-to-maximize-a-plt-show-window-using-python '''
# print '#1 Backend:',plt.get_backend()

# http://stackoverflow.com/questions/12439588/how-to-maximize-a-plt-show-window-using-python
# print '#3 Backend:',plt.get_backend()


colors = ['aliceblue', 'aqua', 'blue', 'chocolate', 'cyan','darkgreen','green',
          'gray','red','ivory','indigo','lavender',
          'lightseagreen', 'olive', 'orange', 'purple','salmon']


if not mbs.stat_files_to_dic():
    exit()
    
mbs.pstat_critic()
    
iPar = 0
numCols = 8
numLines = 16
    
for par in params:
    iPar += 1
    summary[par] = {}
    
    #fig = figList[iPar]
    fig = plt.figure(iPar)
    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.90)
    mng = plt.get_current_fig_manager()

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
                    
    elif plt.get_backend() == 'wxAgg': 
        mng.frame.Maximize(True)
    elif plt.get_backend() == 'QT4Agg': 
        mng.window.showMaximized()
        
    seqXs = []
    mini, maxi = float('inf'), -float('inf')
    
    for min_max in mbs.min_maxs:
        for alig_cons in mbs.alig_conss:
            for std_covar in mbs.std_covars: 
                   
                study = min_max + '-' + alig_cons + '-' + std_covar
             
                try:
                    vals = mbs.dic_detaild_params[study][par]
                    if par == "LnL":
                        vals = [np.log10(-x) for x in vals]

                    
                    if min(vals) < mini:
                        mini = min(vals)
                    elif max(vals) > maxi:
                        maxi = max(vals)                  
                except:
                    pass    
    if par not in ["LnL"]:
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
                            
    print ticks
    
    cntCols = 0
    stri_summary = "   \tAligned\tConsensus\n"
    for min_max in mbs.min_maxs:
        for alig_cons in mbs.alig_conss:
            
            listData = []
            print_ok = False
            
            cntLine = 0
            
            for std_covar in mbs.std_covars: 
                   
                study = min_max + '-' + alig_cons + '-' + std_covar
                
                dic = mbs.dic_detaild_params[study] 
                summary[par][study] = {}
              
                try:
                    summary[par][study] = {}
                    vals = dic[par]
                    N = len(vals)
                    seqXs.append(vals)

                    # mu,muLeft,muRight = mbs.mean_confidence_interval(vals, confidence=0.95)
                    mu, sigma, cum, q025, q1, q2, q3, q975 = mbs.my_ppf([1,2,3,4,5])
                    mu, sigma, cum, q025, q1, q2, q3, q975 = mbs.my_harmonic_ppf([1,2,3,4,5])
                    
                    mu, sigma, cum, q025, q1, q2, q3, q975 = mbs.stat_ppf(vals)
                    

                    summary[par][study]["N"] = N
                    summary[par][study]["mu"] = mu
                    summary[par][study]["sigma"] = sigma
                    summary[par][study]["SE"] = sigma/sqrt(N)
                    summary[par][study]["min"] = min(vals)
                    summary[par][study]["max"] = max(vals)
                    summary[par][study]["median"] = q2
                    summary[par][study]["q025"] = q025
                    summary[par][study]["q975"] = q975

                    try:
                        ''' LnL doesn't have ESS and PSRF - adopted from TL
                            total tree length (the sum of all branch lengths, TL) '''
                        if par == "LnL":
                            summary[par][study]["avgESS"] = mbs.dic_pstat[study]["TL"]['avgESS']
                        else:
                            summary[par][study]["avgESS"] = mbs.dic_pstat[study][par]['avgESS']  
                    except:
                        summary[par][study]["avgESS"] = 0
                    try:
                        if par == "LnL":
                            summary[par][study]["PSRF"] = mbs.dic_pstat[study]["TL"]['PSRF']
                        else:
                            summary[par][study]["PSRF"] = mbs.dic_pstat[study][par]['PSRF']
                    except:
                        summary[par][study]["PSRF"] = 10000
                
                    if summary[par][study]["avgESS"] < 90:
                        sAvgESS = "avgESS=%3.1f **"%summary[par][study]["avgESS"]
                    else:
                        sAvgESS = "avgESS=%3.1f"%summary[par][study]["avgESS"]
                        
                    if (summary[par][study]["PSRF"] - 1) > .1:
                        sPSRF = " PSRF=%1.2f ***"%summary[par][study]["PSRF"]
                    else:
                        sPSRF = " PSRF=%1.2f"%summary[par][study]["PSRF"]

                    if cntLine == 0:
                        ax = plt.subplot2grid((numLines,numCols), (0, cntCols), rowspan=4,colspan=2)
                    else:
                        ax = plt.subplot2grid((numLines,numCols), (6, cntCols), rowspan=4,colspan=2)
                      
                 
                    ''' ax = fig.add_subplot(numLines,numCols,fig_index[numFigure]) '''

                    print_ok = True

                    if par == "LnL":
                        stri_summary = mbs.show_lnl( plt, ax, dic, study, mu, sigma, q2, q025, q975, sAvgESS, sPSRF, vals, cntLine, cntCols, min_max, alig_cons, stri_summary)
        
                    else:
                        mbs.show_distrib(plt, ax, par, study, mu, sigma, q2, q025, q975, sAvgESS, sPSRF, ticks, vals, colors, iPar, mini, maxi, cntCols)
        
   
                    listData.append(vals) 
                    cntLine += 1
                except:
                    pass

            if print_ok:    
                if len(listData) == 2:
                    mbs.show_conf_inteval(min_max, alig_cons, par, summary, plt, numLines,numCols, cntCols, ticks, ticks0, mini, maxi, listData)
                    mbs.show_qqplot2(min_max, alig_cons, par, summary, plt, numLines,numCols, cntCols, listData)
                else:
                    pass
       
            cntCols += 2          
    
    if par == "LnL":
        stri = ""
    else:    
        f_value, p_value = mbs.calc_anova(seqXs)
        if p_value <= 0.05:
            stri = ', at least one distribution is statistically different.'
        else:
            stri =  ', the distributions are statistically similar.'
            
        stri += 'ANOVA: f_value %f   p_value %.5e' %(f_value, p_value)

    left = .05
    top = .97
    fig.text(left, top, par+stri, color="red") 


'''
try:
    #plt.switch_backend('TkAgg') #TkAgg (instead Qt4Agg)
    mng = plt.get_current_fig_manager()
    ### works on Ubuntu??? >> did NOT working on windows
    # mng.resize(*mng.window.maxsize())
    mng.window.state('zoomed') #works fine on Windows!    
except:
    pass
'''
    
print stri_summary
plt.show()
exit()



stri = "%s \t %i \t %2.3e \t %2.3e \t %2.3e \t %2.3e \t %2.3e \t %2.3e \t %s \t %2.3e \t %2.3e"

for par in params:
    print "Param: %s"%(par)
    dic = summary[par] 
    lista = dic.keys()
    lista = sorted(lista)
    
    print "study".ljust(30) + " \t N   \t mu         \t sigma    \t SE          \t minimum  \t maximum  \t median   \t conf.interval 95% \t avgESS \t PSRF"
    

    for min_max in mbs.min_maxs:
        for alig_cons in mbs.alig_conss:
            for std_covar in mbs.std_covars: 
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
         
         

mbs.pstat_critic()

mbs.write_stat_files()
    


'''
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
          
                    



        
        
            
            
    