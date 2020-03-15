import sys
from datetime import datetime

def logStart(inputList, dataList):
   """Creates the log file 'ComplexBuilderLog' and writes initial message in log file
      with the date and the input files"""
   now = datetime.now() 
   dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
   sys.stderr = open("./"+dataList[2]+"/ComplexBuilder.log", "w")
   sys.stderr.write("\n                              ComplexBuilder\n")
   sys.stderr.write("     BUILDING OF COMPLEXES USING THE PAIRED INTERACTION OF ITS ELEMENTS\n")
   sys.stderr.write("Written by Paula Gomis Rosa, Arturo González Vilanova and Marta López Balastegui\n\n")
   print("Job starting time: ", dt_string, file=sys.stderr)
   sys.stderr.write("\nInput files:\n\t - FASTA sequences in file:\n\t ")   
   fas=inputList[-1]
   print(fas, file=sys.stderr)
   sys.stderr.write("\n\t - PDB files:\n\t ")  
   pdbs=list( inputList[i] for i in range(len(inputList)-1) )
   print(pdbs, file=sys.stderr)

def logPairs(inputPairs):
   """Write in the log file the pairing of the complex"""
   return None
