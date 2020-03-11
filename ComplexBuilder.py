import argparser
import logProgress
import interface
from utilities import *

# Read the command line arguments and check format
options = argparser.readArgs()
if options.gui: #enter files through GUI
    inputList=interface.initGui()
    inputFiles = checkInputs(inputList[0], inputList[1])
    if  inputList[2]:
        logProgress.logStart(inputFiles)
else: #enter files through command line
    inputFiles = checkInputs(options.infasta, options.inpdb)
if options.verbose: # if verbose is ON write progress in "ComplexBuilder.log"
    logProgress.logStart(inputFiles)

# Store information of the complex pairs
pdb_files = inputFiles[:-1]
fasta_file = inputFiles[-1]
a = data_extraction(pdb_files,fasta_file)
b = seq_dictionary(a)
print(b)
