import os
from os.path import isfile

def checkInputs(fastaFile, PDBDir):
    """Check if the FASTA file and PDB directory is introduced in 
       the command line"""
    if isfile(fastaFile)==False:
        print ("You should introduce the name of the FASTA file "
               "containing the sequences of the elements of the " 
               "complex you want to build after '-fa'.\n"
               "(type -h for more information of the required "
               "format) ")
    if PDBDir==None or os.path.exists(PDBDir)==False:
        print("You should introduce the name of the directory "
               "containing the PDB files for the element pairs "
               "of the complex you want to build after '-pdb'.\n"
               "(type -h for more information of the required "
               "format) ")
    return None