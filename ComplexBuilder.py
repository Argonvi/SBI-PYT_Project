import argparser
import utilities
import logProgress

# Read the command line arguments and check format
options = argparser.readArgs()
inputFiles = utilities.checkInputs(options.infasta, options.inpdb)
if options.verbose: # if verbose is ON write progress in "ComplexBuilder.log"
    logProgress.logStart(inputFiles)

