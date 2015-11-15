'''
Created on Aug 31, 2015

@author: flavio
'''
import os, glob

ori = ["Drosophila_maxmer_Gene_Adh_consensus_covarion.nxs",\
       "Drosophila_maxmer_Gene_Adh_aligned_covarion.nxs",\
       "Drosophila_mincut_Gene_Adh_consensus_covarion.nxs",\
       "Drosophila_mincut_Gene_Adh_aligned_covarion.nxs"]

L = len("Drosophila_maxmer_Gene_Adh_")

by_these = ["Drosophila_maxmer_Gene_Adh_100L_cutoff7_",\
        "Drosophila_maxmer_Gene_Adh_100L_cutoff7_",\
        "Drosophila_mincut_Gene_Adh_100L_cutoff7_",\
        "Drosophila_mincut_Gene_Adh_100L_cutoff7_"]


path = "/home/flavio/mrbayes/"
path = "c:/mrbayes/"


print "-----------------"

os.chdir(path)
for filename in glob.glob("*.nxs*"):
    for i in range(4):
        if filename.find(ori[i]) == 0:
            print(filename),
            newName = by_these[i] + filename[L:]
            print " found, and replace", filename
            
            os.rename(filename, newName)
            break

    
                        