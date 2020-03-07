import sys
from datetime import datetime

def logStart(inputList):
    """Creates the log file 'ComplexBuilderLog' and writes initial message in log file
       with the date and the input files"""
    now = datetime.now() 
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    sys.stderr = open("./ComplexBuilder.log", "w")
    sys.stderr.write("\n                         ComplexBuilder\n")
    sys.stderr.write("  BUILDING OF COMPLEXES USING THE PAIRED INTERACTION OF ITS ELEMENTS\n")
    sys.stderr.write("Writen by Paula Gomis Rosa, Arturo González and Marta López Balastegui\n\n")
    print("Job starting time: ", dt_string, file=sys.stderr)
    sys.stderr.write("Input files:\n -PDB files: ")
    pdbs=list( inputList[i] for i in range(len(inputList)-1) )
    print(pdbs, file=sys.stderr)
    fas=inputList[-1]
    sys.stderr.write(" -FASTA sequences in file:\n ")
    print(fas, file=sys.stderr)
         