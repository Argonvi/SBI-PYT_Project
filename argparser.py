import sys,os,argparse, gzip
from os.path import isfile, join




def readArgs():
    """ Read and organize the command line arguments and return the namespace"""
    parser = argparse.ArgumentParser(description="Build a protein complex from a set "
                                                 "of pdb files containing the paired "
                                                 "structures of the monomers.")

    parser.add_argument('-fa', '--fasta', dest = "infasta", action = "store", default = None, 
                        help = """FASTA file with the sequences of the proteins
                               conforming the complex you want to build.""")

    parser.add_argument('-pdb', '--pdbDir', dest = "inpdb", action = "store", default = None, 
                        help = """PDB diretory containing the PDB files with the 
                               structure of the pairs conforming the complex you 
                               want to build.""")
                        

    parser.add_argument('-o', '--output', dest = "outfile", action = "store", default = "ComplexBuilding", 
                        help = """Directory where the complex results will be stored. 
                                If it is not defined a new directory 'ComplexBuilding' 
                                will be created.""")

    parser.add_argument('-v', '--verbose', dest = "verbose", action = "store_true", default = False, 
                        help = "Show the detailed progression of the building.")

    parser.add_argument('-st', '--stoichiometry', dest = "stoich", action = "store_true", default = False, 
                        help = "Define a determined stoichiometry.") 

    parser.add_argument('-gui', '--graphicInterface', dest="gui",action="store_true",default=False,
                        help="Display the graphic iterface") 

    return parser.parse_args()




