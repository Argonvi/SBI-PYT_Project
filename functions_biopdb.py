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

#Performing aligment

K = parser.get_structure("K","6gmh_K.pdb")
E = parser.get_structure("E","6gmh_E.pdb")

K_seq = ppb.build_peptides(K)[0].get_sequence()
E_seq = ppb.build_peptides(E)[0].get_sequence()

alignment = pairwise2.align.globalxx(K_seq, E_seq) #The aligment object is a list of the best possible aligments between the sequences
#The default scoring system is (match = +1, mismatch = 0, gap opening and extension = -1). This means that if the score is the same as the sequence length, both sequences are identical.
print(pairwise2.format_alignment(*alignment[0])) #We use a special function to print the results

alignment_s = pairwise2.align.globalxx(K_seq, E_seq, score_only = True) #This returns only the scores
