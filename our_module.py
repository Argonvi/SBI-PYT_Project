from Bio.PDB import *
from Bio import SeqIO, pairwise2


def inter(pdb_files, fasta_file, threshold = 0.95):
    """Given a list of pdb file paths and a fasta file path (with the sequences in the full complex)
    , returns a dictionary with fasta sequence ids as keys, and a tuples as values. The tuple's
    first value is the file where the sequence has been found and the second value is the chain id
    in such file.
    The threshold of similarity considered for the sequences to be the same is 0.95 at default."""
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
                    inter[fasta_id].append((pdbfile.split("/")[-1], chain.get_id()))
    return inter
