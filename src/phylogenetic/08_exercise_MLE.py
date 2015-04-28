'''
Created on 04/02/2015

@author: http://biopython.org/wiki/PAML

'''
from Bio import AlignIO, Phylo
from Bio.Phylo.PAML import codeml

aln = AlignIO.read('msa.phy', 'phylip')
print aln

'''
    codeml : for codons or aa
    baseml : for nucleotides
'''

cml = codeml.Codeml(alignment = "msa.phy", tree = "msa_nj.nhx", out_file = "results.out", working_dir = "./scratch")
''' idem
cml = codeml.Codeml()
cml.alignment = "align.phylip"
cml.tree = "species.tree"
cml.out_file = "results.out"
cml.working_dir = "./scratch"

'''

'''
 Options may be set by the set_option() function and 
 their values may be retrieved by the get_option() function:
'''


cml.set_options(clock=1)
cml.set_options(seqtype=1)  # 1 codon seq ~ exon, 2 for amino acids
cml.set_options(NSsites=[0,1,2])
cml.set_options(aaRatefile="wag.dat")
cml.get_option("NSsites")

opt = False

if opt:
    cml.read_ctl_file("codeml.ctl")
    cml.print_options()
    
'''
Running the program
Executing the object's run() method will run codeml with the current options and will return the parsed results in a dictionary object (see below). You can also pass a number of optional arguments to the run() method:
verbose (boolean): codeml's screen output is suppresed by default; set this argument to True to see all of the output as it is generated. This is useful if codeml is failing and you need to see its error messages.
parse (boolean): set to False to skip parsing the results. run() will instead return None
ctl_file (string): provide a path to an existing control file to execute. The file is not parsed and read into Codeml's options dictionary. If set to None (default), the options dictionary is written to a control file, which is then used by codeml.
command (string): provide a path to the codeml executable. This is set to "codeml" by default, so if the program is in your system path, that should suffice. If it's not in your system path or, for example, if you use multiple versions of PAML, you may instead provide the full path to the executable.
'''

cml.run(ctl_file="codeml.ctl", verbose=True, command="codeml")

# results = codeml.read()

