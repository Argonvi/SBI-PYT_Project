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

ascii_letters_generator = (a for a in string.ascii_letters)

def checkCommands(commands):
    """Check the mode of operation of ComplexBuilder: if '-gui' has been defined
        it must be a unique argument, otherwise use the rest of commandline arguments"""
    inputs=[]
    data=[]
    if commands.gui: #enter files through GUI
        if len(sys.argv)==2:
            inputList=interface.initGui()
            inputs = checkInputs(inputList[0], inputList[1]) #FASTA file, PDB dir
            data.append(inputs) # list of inputs : FASTA + PDB dir, data[0]
            if  inputList[2]: # verbose T/F
                logProgress.logStart(inputs)
            data.append(inputList[3])  # stoichiometry None/Filename, data [1]
            data.append(inputList[4]) #output file data[2]
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
        checkSt(commands.stoich,data)
        checkOutput(commands.outfile, data)

    return data

def checkSt(stName, inputsListed):
    """Check the stoichiometry file existence entered by the user and
    append it to the list of input values, otherwise add None"""
    if stName:
        if isfile(stName) and stName.endswith('.txt'):
            inputsListed.append(stName)
        else:
            raise ValueError("""
                You should introduce the name of an existing .txt
                file with the sctoichiometry data after '-st'.\n
                Type -h for more information of the required
                format""")
    else:
        inputsListed.append(None)
    return None


def checkOutput(outputName, inputsListed):
    """Check the output name existence defined by the user and
    append it to the list of input values, otherwise raise an error"""
    if outputName is not None:
        inputsListed.append(outputName)
    else:
        raise ValueError("""
                You should introduce the name of the directory
                where the resulting files will be stored after
                 '-o'.\n
                Type -h for more information of the required
                format""")
    return None

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
        if len(inputList)<1: # not PDB files inside PDBdir
            raise ValueError("""
                You should introduce the name of the directory
                containing the PDB files for each interacting
                pair of the complex after '-pdb'.\n
                Type -h for more information of the required
                format""")
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
    interactions wich have:
     1: The model in which the sequence is found
     2: The corresponding chain in said model
     3: The other sequence (id) with which it is interacting"""
    sequence_dictionary = {}
    for k, v in data.items():
        for k2, v2 in v.items():
            sequence_dictionary.setdefault(v2, [])
            if list(k.get_chains()) in [list(a[0].get_chains()) for a in sequence_dictionary[v2]]: continue
            if list(v.values())[0]==list(v.values())[1]: other_seq=list(v.values())[0]
            else: other_seq = [seq_id for seq_id in v.values() if seq_id != v2][0]
            sequence_dictionary[v2].append((k,k[k2],other_seq))
    return sequence_dictionary

def stoichometry(file, information):
    """Takes as input a file with the stoichomety information and stores
    it in a dictionary format with sequences as keys and number of appearances
    in the complex as values.
    The sequence names should be the same as in the fasta file."""
    fasta_ids = list(information.keys())
    dictionary = {}
    if file is not None:
        with open(file) as f:
            for line in f:
                line = line.strip()
                line = line.split(":")
                dictionary[line[0]] = int(line[1])
    for fasta_id in fasta_ids:
        if fasta_id not in list(dictionary.keys()):
            dictionary[fasta_id] = 1
    return dictionary

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

    char = next(ascii_letters_generator)
    third_chain.id = char
    macrocomplex.add(third_chain)
    #if third chain collapses with another chain, raise an error
    return macrocomplex

def writte_pdb(structure,directory,name_pdb):
    io = PDBIO()
    io.set_structure(structure)
    save=directory+"/"+name_pdb
    io.save(save)


def constructor(information,stoich):
    """Description"""
    chains_in_complex={}

    #Choose a model to start iterating
    for seq, interactions in information.items():
        if len(interactions)>1 and (interactions[0][2]!=interactions[1][2] or stoich[interactions[0][2]]>1): break
    rand_interaction = random.choice(information[seq])
    start_model = rand_interaction[0]
    first_chain = rand_interaction[1]

    #Store its chains in the complex dictionary and used chains list
    chains_in_complex[seq] = [first_chain]
    chains_in_complex[rand_interaction[2]] = [chain for chain in start_model.get_chains() if chain.get_id() != first_chain.get_id()]
    chains_used = [seq]

    #Make a copy to modify
    start_model_copy = copy.deepcopy(start_model)

    #Extend this model as much as possible from the first chain
    for interaction in information[seq]:
        if interaction[0] is start_model: continue
        second_model = interaction[0]
        same_chain = interaction[1]
        other_seq = interaction[2]
        third_chain = [chain for chain in second_model.get_chains() if chain.get_id() != same_chain.get_id()][0]
        complex_out=superimpositor(first_chain, same_chain, third_chain, start_model_copy)
        chains_in_complex.setdefault(other_seq,[])
        chains_in_complex[other_seq].append(third_chain)

    #From the resulting model, keep adding chains, until the model has as many chains as specified in the stoichiometry
    while len(list(complex_out.get_chains())) < sum(stoich.values()):
        #Get a sequence that is in the complex but is yet to be used as a core for the extension of the complex
        seq = [chain for chain in chains_in_complex if chain not in chains_used][0]

        for first_chain in chains_in_complex[seq]:
            chains_used.append(seq)

            for interaction in information[seq]:
                if (interaction[2] in chains_in_complex) and (len(chains_in_complex[interaction[2]])==stoich[interaction[2]]): continue
                second_model = interaction[0]
                same_chain = interaction[1]

                for chain in second_model.get_chains():
                    if chain.get_id() == same_chain.get_id(): continue
                    complex_out=superimpositor(first_chain, same_chain, chain, complex_out)
                    chains_in_complex.setdefault(interaction[2],[chain])
                    chains_in_complex[interaction[2]].append(chain)
    return complex_out
