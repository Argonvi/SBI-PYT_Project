import os
from os.path import isfile, join

from Bio.PDB import *
from Bio import SeqIO, pairwise2
import random

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

def data_extraction(pdb_files, fasta_file, threshold = 0.95):
    """Takes as input a list of pdb files and a fasta file
    Returns a dictionary of dictionaries. The primary key is the model, the secondary
    key is the chain_id in said model and the value is the fasta_id of said chain."""
    big_dictionary = {}
    for pdb_file in pdb_files:
        model = PDBParser().get_structure(pdb_file.split(".")[0],pdb_file)[0]
        big_dictionary[model] = {}
        for chain in model.get_chains():
            pdb_seqs = PPBuilder().build_peptides(chain)
            if len(pdb_seqs)>1:
                pp_seq = "".join(list([str(pp.get_sequence()) for pp in pdb_seqs]))
            else:
                pp_seq = pdb_seqs[0].get_sequence()
            for seq_record in SeqIO.parse(fasta_file,"fasta"):
                fasta_seq = seq_record.seq
                fasta_id = seq_record.id
                score = pairwise2.align.globalxx(fasta_seq,pp_seq, score_only = True)
                normalized_score = score/len(max([fasta_seq,pp_seq]))
                if normalized_score > threshold:
                    big_dictionary[model][chain.get_id()] = fasta_id
                    continue
    return big_dictionary

def seq_dictionary(data):
    """Transforms the dictionary of dictionaries given by data_extraction()
    into a dictionary with the fasta_ids as keys and as values, a list of
    tuples wich have:
     1: The model in which the sequence is found
     2: The corresponding chain in said model
     3: The other sequence (id) with which it is interacting"""
    sequence_dictionary = {}
    for k, v in data.items():
        for k2, v2 in v.items():
            sequence_dictionary.setdefault(v2, [])
            for seq_id in v.values():
                if seq_id != v2: other_seq = seq_id
            sequence_dictionary[v2].append((k,k[k2],other_seq))
    return sequence_dictionary

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
