import os
from os.path import isfile, join

from Bio.PDB import *
from Bio import SeqIO, pairwise2
import random
import copy
 

def checkInputs(fastaFile, PDBDir):
    """Check if the FASTA file and PDB directory is introduced in
       the command line and returns all the files in a list"""

    if isfile(fastaFile)==False:
        raise ValueError("""You should introduce the name of the FASTA file
               containing the sequences of the elements of the
               complex you want to build after '-fa'.\n
               (type -h for more information of the required
               format) """)
    if PDBDir==None or os.path.exists(PDBDir)==False:
        raise ValueError("""You should introduce the name of the directory
               containing the PDB files for the element pairs
               of the complex you want to build after '-pdb'.\n
               (type -h for more information of the required
               format) """)
    if isfile(fastaFile) and os.path.exists(PDBDir):
        if PDBDir.endswith("/"):
            string = ""
        else:
            string = "/"
        inputList=[]
        inputList = [PDBDir+string+f for f in os.listdir(PDBDir) if f.endswith(".pdb")  and isfile(join(PDBDir, f))]
        inputList.append(fastaFile)
        return inputList
    else:
        return None

def seq_finder(pdb_files, fasta_file, threshold = 0.95):
    """Given a list of pdb file paths and a fasta file path (with the sequences in the full complex)
    , returns a dictionary which has has fasta sequence ids as keys, and a tuples as values. The tuple's
    first value is the model to which the chain belongs and the second value is the chain object
    in such model."""
    inter = {}
    for seq_record in SeqIO.parse(fasta_file,"fasta"):
        fasta_seq = seq_record.seq
        fasta_id = seq_record.id
        inter[fasta_id] = []
        for pdbfile in pdb_files:
            pdb_data = PDBParser().get_structure(pdbfile.split(".")[0],pdbfile)[0]
            for chain in pdb_data:
                pdb_seqs = PPBuilder().build_peptides(chain)
                if len(pdb_seqs)>1:
                    pp_seq = "".join(list([str(pp.get_sequence()) for pp in pdb_seqs]))
                else:
                    pp_seq = pdb_seqs[0].get_sequence()
                score = pairwise2.align.globalxx(fasta_seq,pp_seq, score_only = True)
                normalized_score = score/len(max([fasta_seq,pp_seq]))
                if normalized_score > threshold:
                    inter[fasta_id].append((pdb_data, chain))
    return inter


def superimpositor(first_chain, same_chain, third_chain,macrocomplex):
    """REVISAR ESTA DESCRIPCION  PORQUE ES UN CHURRO"""
    
    """This function takes as input 3 chain objects:
            The chain that we take as reference in order to do the superimposition
            The chain that we want to superimpose 
            The chain that we are going to rotate in because it interacts with the previous chain
        It also takes the macrocomplex in order to add the third chain. 
    """ 
    atom_list1 = Selection.unfold_entities(first_chain, 'A')
    atom_list2 = Selection.unfold_entities(same_chain, 'A')
    sup = Superimposer()
    sup.set_atoms(atom_list1, atom_list2)
    sup.apply(third_chain)
    macrocomplex.add(third_chain)
    return macrocomplex




def constructor(information):
    #Get a core model randomly
    
    #desde aqui 
    rand_seq = random.choice(list(information.keys()))
    rand_tupla = random.choice(information[rand_seq])
    rand_model=rand_tupla[0]
    first_chain =rand_tupla[1]
    rand_model2 = copy.deepcopy(rand_model)
    
    for tupla in information[rand_seq]:
        if tupla[0] is rand_model:
            continue
        second_model = tupla[0]
        same_chain = tupla[1]
        for chain in second_model.get_chains():
            if chain.get_id() != same_chain.get_id():
                third_chain = chain
        rand_model2=superimpositor(first_chain, same_chain, third_chain,rand_model2)
        
    return rand_model2    


    #hasta aqui hemos conseguido unas seqs random 
##first round           
    
    
  
        
        
        
