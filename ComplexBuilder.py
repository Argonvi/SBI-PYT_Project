import argparser
import logProgress
import interface
from utilities import *


# Read the command line arguments and check format
options = argparser.readArgs()
inputData = checkCommands(options)  


# Store information of the complex pairs
inputFiles=inputData[0]
pdb_files = inputFiles[:-1]
fasta_file = inputFiles[-1]
st_file=inputData[1] # if there is not stoichiometry = None
out_file=inputData[2] 

a = data_extraction(pdb_files,fasta_file)
b = seq_dictionary(a)
#print(b)

mycomplex=constructor(b,stoich)

writte_pdb(mycomplex,out_directory,"file_name.pdb")
