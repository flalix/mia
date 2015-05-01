---
title: "MIA"
author: "Flavio Lichtenstein"
date: "Wednesday, April 29, 2015"
output: html_document
---

## How to install

There are two ways to install:  

  - clone the project via github  
  - download the project 
  

## How to clone?

Here are the commands to clone MIA project:  

  - Install GitHub application  
    + for windows: [click here](https://github.com/)  
    + for linux give the following command: **sudo apt-get install git** 

<br />

  - Run GitHub and write the following commands:
  
git init
git clone https://github.com/flalix/mia/



## How to download? 

  - Go to [MIA site](https://github.com/flalix/mia/) and click in download
  - Click in **download** 


![MIA site in GitHub](https://github.com/flalix/mia/blob/master/image/github_mia_site.png?raw=true)

<br />


  - Create a directory in windows:

cd c:\\Users\\your_name  
md mia  

<br />


  - Create a directory in linux:
    
cd ~  
mkdir mia  

<br />

  - and **decompress the zip file here in mia directory**

<br />

## Running MIA

Once cloned or dowloaded (and unziped) MIA archives, you will find the executable file (binary) at:  

  - if windows:  c:\\user\\your_name\\mia\\exe\\windows  
  
  - if linux: ~/mia/exe/linux 
  
Give a double click in MIA and start it.  


<br />



## First time at front end and what to do?

If you stared MIA correctly this front end will appear:

![MIA front end](https://github.com/flalix/mia/blob/master/image/mia_first_time.png?raw=true)

<br />

Now you have four alternatives to start your analysis:  

1. Click in [Get GBK from NCBI] and you will download all gbk from the organism "Drosophila" and gene "Adh".  
2.	Change "Adh" to "AMY" and in Gene List write "AMY, AMYREL", and click in [Get GBK from NCBI] and you will analyze a short and interesting sample data.  
3.	You also may define your own Organism and Gene. If necessary you may define words in Title to find any desired gene-experiment.
4.	Another possible way is to see a previous analysis that we did. In mia/data you will find 2 set of samples. One for "Adh" and the other for "AMY".  

  - click in "save" and exit the MIA program  
  - look for:
    + in windows:  c:\\user\\your_name\\Drosophila\\data
    + in linux: ~/Drosophila/data  
    
  - you will find /data/Adh  
  
  - create also /data/AMY  
  - go to your clone/download directory and find data/adh and data/AMY  
  
  - copy each fasta.rar and files.rar to its respective directory  
  
  - decompress and you will find:  
    + ../ Drosophila/data/Adh/fasta  
    + ../ Drosophila/data/Adh/files  
    
  - both directories have many files inside  
  
  - in mia/data you will find **"default.ini"**, **overwrite** it in ../ Drosophila/ data  
  
  - find again MIA executable (binary) and **start it**  
  
  - in the first combo box you will find " Drosophila - Adh" and " Drosophila - AMY".  
  - choose one of these options  
  - jump (click) in VMI, or HMI, or JSD or Cluster tab and see the results  
  - any other information sarch in the MIA manual, please  
  







