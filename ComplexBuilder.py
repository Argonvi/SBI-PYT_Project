import argparser
from utilities import *
import logProgress


# Read the command line arguments and check format
options = argparser.readArgs()
inputFiles = checkInputs(options.infasta, options.inpdb)
if options.verbose: # if verbose is ON write progress in "ComplexBuilder.log"
    logProgress.logStart(inputFiles)

pdb_files = inputFiles[:-1]
fasta_file = inputFiles[-1]
a = seq_finder(pdb_files,fasta_file)

print(a)

constructor(a)
