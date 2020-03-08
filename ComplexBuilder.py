import argparser
import utilities
import logProgress
from our_module import *

# Read the command line arguments and check format
options = argparser.readArgs()
inputFiles = utilities.checkInputs(options.infasta, options.inpdb)
if options.verbose: # if verbose is ON write progress in "ComplexBuilder.log"
    logProgress.logStart(inputFiles)

pdb_files = inputFiles[:-1]
fasta_file = inputFiles[-1]
print(inter(pdb_files,fasta_file))
