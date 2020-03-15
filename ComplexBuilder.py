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
out_directory=inputData[2] 

a = data_extraction(pdb_files,fasta_file)
b = seq_dictionary(a)
print(b)


stoich=stoichometry(st_file,b)
mycomplex=constructor(b,stoich)
list(mycomplex.get_chains())

writte_pdb(mycomplex,out_directory,"final_prueba.pdb")
 



