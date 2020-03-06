import argparser
import utilities

# Read the command line arguments and check format
options = argparser.readArgs()
utilities.checkInputs(options.infasta, options.inpdb)

