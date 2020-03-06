#This file is for the record of useful functions for the project that can be found in the Bio.PDB package
# link to the documentation: https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ

from Bio.PDB import *

#Get the pdb structure into a python object
parser = PDBParser()
structure = parser.get_structure("protein_complex","6gmh.pdb")

#Get the polypeptides

ppb = PPBuilder()
polypeptides = ppb.build_peptides(structure) #is a list of polypeptide objects

#Get the sequences

for pp in polypeptides:
    pp.get_sequence() #is a sequence object, similar to a string
