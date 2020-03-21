import argparser
import logProgress
import interface
from utilities import *


# Read the command line arguments and check format
options = argparser.readArgs()
inputData = checkCommands(options)  


# Store input information 
inputFiles=inputData[0]
pdb_files = inputFiles[:-1]
fasta_file = inputFiles[-1]
st_file=inputData[1] # if there is not stoichiometry = None
out_directory=inputData[2] 

# Identify common elements between different pairs
a = data_extraction(pdb_files,inputData[3],fasta_file)
b = seq_dictionary(a)
print(b)

# Add possible stoichimetry 
stoich=stoichometry(st_file,b)
if inputData[3]:
    logProgress.progress(out_directory,b,stoich)

# Build the structure
mycomplex=constructor(b,stoich,inputData[3])
list(mycomplex.get_chains())


# Write results in output directory
path = out_directory+"/resultado.pdb"
write_pdb(mycomplex,path)
if inputData[3]:
    logProgress.end(out_directory)
