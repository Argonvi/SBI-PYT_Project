import os
from os.path import isfile, join
import sys
from Bio.PDB import *
from Bio import SeqIO, pairwise2
import random
import copy
import string
import interface
import logProgress

def checkCommands(commands):
    """Check the mode of operation of ComplexBuilder: if '-gui' has been defined
        it must be a unique argument, otherwise use the rest of commandline arguments"""
    inputs=[]
    data=[]
    if commands.gui: #enter files through GUI
        if len(sys.argv)==2:
            inputList=interface.initGui()
            inputs = checkInputs(inputList[0], inputList[1]) #FASTA file, PDB dir
            data.append(inputs)
            if  inputList[2]: # -v selected in GUI
                logProgress.logStart(inputs)
            if inputList[3]:
                data.append(inputList[3])
        else:
            raise ValueError("""
                If you want to use the graphical interface to 
                introduce the complex data, the only commandline 
                argument possible is '-gui'.\n
                Type -h for more information of the required
                format""")
    else: #enter files through command line
        inputs = checkInputs(commands.infasta, commands.inpdb)
        data.append(inputs)
    if commands.verbose: # if verbose is ON write progress in "ComplexBuilder.log"
        logProgress.logStart(inputs)
    if commands.stoich:
        data.append(commands.stoich)
    else:
        data.append(None)
    
    print('utilities inputs',inputs)
    print('utilities data',data)
    return data

def checkInputs(fastaFile, PDBDir):
    """Check if the FASTA file and PDB directory is introduced in
       the command line and returns all the files in a list"""

    if isfile(fastaFile)==False:
        raise ValueError("""
                You should introduce the name of the FASTA file
                containing the sequences of the elements of the
                complex you want to build after '-fa'.\n
                Type -h for more information of the required
                format""")
    if PDBDir==None or os.path.exists(PDBDir)==False:
        raise ValueError("""
                You should introduce the name of the directory
                containing the PDB files for each interacting
                pair of the complex after '-pdb'.\n
                Type -h for more information of the required
                format""")
    if isfile(fastaFile) and os.path.exists(PDBDir):
        if PDBDir.endswith("/"):
            string = ""
        else:
            string = "/"
        inputList=[]
        inputList = [PDBDir+string+f for f in os.listdir(PDBDir) if f.endswith(".pdb") and isfile(join(PDBDir, f))]
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
            if list(v.values())[0]==list(v.values())[1]:
                other_seq=list(v.values())[0]
            else:
                for seq_id in v.values():
                    if seq_id != v2: other_seq = seq_id
            sequence_dictionary[v2].append((k,k[k2],other_seq))
    return sequence_dictionary



def stoichometry(file, information):
    """Takes as input a file with the stoichomety information and stores
    it in a dictionary format with sequences as keys and number of appearances
    in the complex as values.
    The sequence names should be the same as in the fasta file."""
    fasta_ids = list(information.keys())
    dictionary = {}
    with open(file) as f:
        for line in f:
            line = line.strip()
            line = line.split(":")
            dictionary[line[0]] = line[1]
    for fasta_id in fasta_ids:
        if fasta_id not in list(dictionary.keys()):
            dictionary[fasta_id] = 1
    return dictionary
        


def len_complex(stoich):
    length=0
    for i in stoich:
        length+=stoich[i]
    return length

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
    print(first_chain.get_id(), same_chain.get_id())
    print(len(atom_list1))
    print(len(atom_list2))

    sup.set_atoms(atom_list1, atom_list2)
    sup.apply(third_chain)
    try:
        macrocomplex.add(third_chain)
    except:
        third_chain.id=random.choice(string.ascii_letters)
        macrocomplex.add(third_chain)
    return macrocomplex


def constructor(information,stoich):
    #Get a core model randomly
    chains_in_complex={}
    chains_dict_used=[]
    n=0
    
    #first 
    for seq in information:
        if len(information[seq])>1 and (information[seq][0][2]!=information[seq][1][2] or stoich[information[seq][0][2]]>1 ):
            rand_seq=seq
            break
    rand_tupla = random.choice(information[rand_seq])
    rand_model=rand_tupla[0]
    first_chain =rand_tupla[1]
    chains_in_complex[rand_seq]= [first_chain]
    n+=1
    chains_dict_used.append(rand_seq)
    rand_model2 = copy.deepcopy(rand_model)
    for chain in rand_model.get_chains():
        if chain.get_id() != first_chain.get_id():
            chains_in_complex[rand_tupla[2]]=[chain]
            n+=1
    
        
    for tupla in information[rand_seq]:
        if tupla[0] is rand_model:
            continue
        second_model = tupla[0]
        same_chain = tupla[1]
        
        for chain in second_model.get_chains():
            if chain.get_id() != same_chain.get_id():
                third_chain = chain
                rand_model2=superimpositor(first_chain, same_chain, third_chain,rand_model2)
                chains_in_complex[tupla[2]]=[third_chain]
                n+=1
    
    ##following  
    
    while n<len_complex(stoich):

        for chain in chains_in_complex:
            if chain not in chains_dict_used:
                rand_seq=chain
                break
    
        for i in range(len(chains_in_complex[rand_seq])):
            first_chain =chains_in_complex[rand_seq][i-1]
            chains_dict_used.append(rand_seq)
    
            for tupla in information[rand_seq]:
                if ((tupla[2] in chains_in_complex) and (len(chains_in_complex[tupla[2]])==stoich[tupla[2]])):
                    continue
                second_model = tupla[0]
                same_chain = tupla[1]
            
                for chain in second_model.get_chains():
                    if chain.get_id() != same_chain.get_id():
                        third_chain = chain
                        rand_model2=superimpositor(first_chain, same_chain, third_chain,rand_model2)
                        if tupla[2] in chains_in_complex:
                            chains_in_complex[tupla[2]].append(third_chain)
                            n+=1
                        else:
                            chains_in_complex[tupla[2]]=[third_chain]
                            n+=1
   
    
    return rand_model2   
