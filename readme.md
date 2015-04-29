---
title: "MIA"
author: "Flavio Lichtenstein"
date: "Wednesday, April 29, 2015"
output: html_document
---

### Mutual Information Analyzer, a software that clusters molecular sequences based on entropy and mutual information
  
Federal University of Sao Paulo (UNIFESP), DIS-Bioinformatics    
  
**Motivation:** Here we propose a method to discriminate closely related species at the molecular level using entropy and mutual information. Sequences of orthologous genes in the same species might contain species-specific covariation patterns that can be identified through mutual information.    

**Summary:** Mutual Information Analyzer (MIA) is a pipeline written in Python with the intent to calculate Vertical Entropy, Vertical and Horizontal Mutual Information. From VH, VMI and HMI distributions, Jensen-Shannon Divergence (JSD) is calculated to estimate the distances between species sequences. Each pair of mutual information distribution distances with their respective standard errors are calculated and stored in distance matrices. These distances between distributions can be presented as histograms or hierarchical cluster dendrograms.    
<br />

#### Methods:
    
**Mutual Information Analyzer (MIA)** is a pipeline written in Pytho with the following algorithms:

- A1)  NCBI: gathers data in NCBI and stores them in GBK  file format; 
- A2) Gbk to Fasta: analyze GBK file and organizes in fasta files per species; 
- A3) Alignment: aligns all sequences and at the end creates two fasta files: "mincut" cutting out columns and sequences with large gaps, and "maxmer" maintaining the maximum possible gaps; 
- A4) Purging: replaces ambiguous nucleotides via IUPAC nucleotide ambiguity table, and eliminates sequences with undesirable words in their names like "synthetic"; 
- A5) Consensus: replaces gaps by their vertical consensus nucleotide; 
- A6) VMI: calculates and stores Vertical Entropy (VH) and Vertical Mutual Information (VMI) distributions, and plots the respective histograms and heat maps; 
- A7) HMI: calculates and stores Horizontal Mutual Information (HMI) distributions, and plots the histograms;
- A8) JSD: calculates Jensen-Shannon Divergence, storing distances and their SEs in distance matrix files, and plots the histograms; 
- A9) HC: calculates hierarchical cluster and present it as a dendrogram; A10) Entropy: simulates Shannon Entropy.
<br />
<br />

HMI and VMI are calculated with and without bias corrections, therefore, the gain or loss of information for "mincut" versus "maxmer", with or without bias correction, can be compared. Distances between distributions are calculated via the square root of JSD. Since Mutual Information and JSD are not linear functions of the data their standard errors are calculated by empirical propagation.

<br />
<br />

####Images:

<br />

**Vertical Shannon Entropy**


![Vertical Entropy](https://github.com/flalix/mia/blob/master/image/VHShannon_Drosophila_mincut_ananassae_Gene_Adh_NOL1_100L_cutoff7.png?raw=true)


<br />

**Vertical Mutual Information 2D Heatmap**

![VMI 2D Heatmap](https://github.com/flalix/mia/blob/master/image/HeatMap_2D_VMI_Drosophila_maxmer_ananassae_Gene_Adh_NOL1_100L_cutoff7_bias_corr.png?raw=true)


<br />



**Vertical Mutual Information 3D Heatmap**

![VMI 3D Heatmap](https://github.com/flalix/mia/blob/master/image/HeatMap_3D_VMI_Drosophila_maxmer_paulistorum_Gene_Adh_NOL1_100L_cutoff7_bias_corr.png?raw=true) 


<br />

  
**Horizontal Entropy**

![HMI](https://github.com/flalix/mia/blob/master/image/HMI_Drosophila_maxmer_ananassae_Gene_Adh_frame0_NOL1_100L_cutoff7.png?raw=true)




<br />

  
**JSD Histogram**

![JSD Histogram](https://github.com/flalix/mia/blob/master/image/JSD_HMI_Drosophila_maxmer_Gene_Adh_frame0_NOL1_100L_cutoff7_bias_corr.png?raw=true)



<br />

  
**Hierarchical Cluster**

![Hierarchical Cluster](https://github.com/flalix/mia/blob/master/image/Cluster_JSD_HMI_Drosophila_maxmer_Gene_Adh_frame0_NOL1_100L_cutoff7_bias_corr.png?raw=true)


