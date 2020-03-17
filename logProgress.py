import sys
import os
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
   return None

def progress(outDir,dic,stValues):
   """Write in the log file the progression information"""
   sys.stderr = open("./"+outDir+"/ComplexBuilder.log", "a")

   keyList=[]
   for element in dic:
      keyList.append(element)

   print("\nNumber of chains found: ", len(keyList), file=sys.stderr) 
   print("IDs of the chains: ", (', '.join(keyList)), file=sys.stderr)

   if stValues is not None:
      print("\nThe stoichiometry for those chains is: ", file=sys.stderr)
      for element in stValues:
         print(element,": ",stValues[element], file=sys.stderr)
   print("\n", file=sys.stderr)
   return None

def clash(boolean,chainName,verboseOn):
   """Show message in logFile when a chain is added or discarted due to clashes"""
   if verboseOn:
      if boolean:
         print("Chain %s has not been added, due to clashes with the macrocomplex structure." %chainName, file=sys.stderr)
      else:
         print("Chain %s has been correctly added." %chainName, file=sys.stderr)
   return None

def end(outDir):
   """ Final message in logFile """
   cwd = os.getcwd()
   resultDir=cwd+'/'+outDir
   print("\nProcess completed correctly",file=sys.stderr)
   print("\nThe results are stored in directory '%s'." %resultDir, file=sys.stderr)
   return None
