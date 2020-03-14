<img src="https://cdn.pixabay.com/photo/2017/10/25/06/13/protein-icon-2887050_960_720.png" title="ComplexBuilderLogo" alt="ComplexBuilderLogo" height="100" width="100">

<!-- [![FVCproductions](https://avatars1.githubusercontent.com/u/4284691?v=3&s=200)](http://fvcproductions.com) -->

<!--***INSERT GRAPHIC HERE (include hyperlink in image)***-->

# Complex Builder

> Generate macrocomplexes superimposing paired interacting elements

> protein - protein or protein - DNA

> the result is stored in a PDB file

## Table of Content

- [Why building complexes?](#whybuildingcomplexes?)
- [Installation](#installation)
- [Options](#options)
- [Examples](#exmples)
- [Limitations](#limitations)
- [Team](#team)


## Why building complexes?

The total number of proteins in humans is around 20K. From this quantity we already know the structure of 4K proteins. However, for 6K we have good templates and for another 6K of them, we have reasonably good templates; that means we have more or less 50-70% of the structure of the human proteins covered.

On the other hand, the number of interactions that may occur between all the proteins is still unknown. There are some studies estimating 130-650K interactions while others estimate more than millions. In this context, 120K interactions are confirmed from experiments and only 7K structures of complexes are known, which is a very small proportion of the total number of possible interactions anyway.

So, the complex-structure coverage, including homologs, is around the 30%, much less than the mentioned coverage of monomers. That means there is still a long way to go into the complex-structure study field. 

Identifying the structure of interacting proteins, complexes, is not an easy task. In monomers, homology modelling can be used to nearly cover all folds to determine a protein structure but, for interactions it is not enough. 

ComplexBuilder tries to generate macrocomplex structures. To do so, the superimposing technique is used: it receives a list of PDB files, each of these files contains the structure of an interacting pair and, by superimposing the common elements of different pairs, it builds the final structure. Although there are other methods to generate macrocomplexes, the superimposition is fast and effective. Furthermore, not only protein-protein interacting pairs can be analyzed, but also proteins with DNA, to end up generating a macrocomplex structure of proteins and DNA chains. 

The core of ComplexBuilder is the `constructor` function, which is the responsible of the building process. The PDBs of the interacting pairs are analyzed and classified using the information of the FASTA file. The common chains in different PDB files are identified: they will have more than 95% of identity in their sequences and an RMSD of less than 7<!-- Confirm this --> angstroms after superimposing their structures. After this classification, `constructor` begins to append elements to the macrocomplex: it starts with a pair, then it looks for another pair with a common element with the starting one and, if after superimposing the common element, the structure has not clashed, it appends the new element. Like this, the macrocomplex have now three elements and, the next step will repeat the process looking for another pair, superimposing the common element and checking the possible clashes. 

<!-- DIAGRAM of the program  --> 

In addition, as a lot of complexes follow a determined stoichiometry, ComplexBuilder can also make the construction following it. If a stoichiometry is defined, the final structure will contain the exact number of elements indicated by it.

Finally, a new PDB file is stored in an output folder(...)

<!-- Chimera comparison with prediction etc  --> 


## Installation

Prerequisites:

- Python 3.0 `https://www.python.org/download/releases/3.0/`
<!-- modeler/itasser/chimera? -->


### Clone

- Clone this repo to your local machine using `https://github.com/Argonvi/SBI-PYT_Project.git`

### Setup

<!-- Needed or it will be already in the package?? --> 
- Install and upgrade BioPython `https://biopython.org/wiki/Download` 

> Using Python package management tool `pip` is easy

```shell
$ pip install biopython
$ pip install biopython --upgrade
```

## Options

ComplexBuilder can be run using command-line arguments or using the graphical interface.



### Command-line

You can introduce the arguments via command-line:

#### Mandatory arguments

To execute ComplexBuilder three arguments are required:

- -fa: FASTA file with the sequences of the proteins or DNA that will conform the complex.

- -pdb: directory containing the PDB files with the structure of the pairs that will conform the complex.

- -o: directory name where the complex results will be stored. 

#### Optional arguments

- -v: show the detailed progression of the building process in a file called 'ComplexBuilder.log'.

- -st: File containing a determined stoichiometry to the complex. The information of the stoichiometry must be: the ID of the sequence chain (concordant with the FASTA file ID) followed by the number of times it has to be present in the complex after ' : '

ID_as_FASTA_file : stoichiometry (one per line) in format .txt. 
Example for a stoichiometry 2A2B for '1GZX':

```shell
1GZXA:2
1GZXB:2
```

### Graphical interface

Otherwise, the macrocomplex can also be build using the graphical interface:

```shell
$ python3 ComplexBuilder.py -gui
```

>In this case just the -gui tag is needed!

- To build the macrocomplex fill in the main window requirements. As for running via command-line, a FASTA file with the sequences and a directory with the PDB files of the interacting pairs are required. In addition, a name for the folder where the results will be stored is needed, after typing it you should confirm it. 

> Note that, as seeing in the following demonstration, to select the PDB directory you have to enter in the desired directory and then select it. 

Furthermore, additional options can be set:

- In the main window you can specify if you want to create a log file where the process of the execution will be displayed. 

- In the top menu, in the 'Options' dropdown, there is the 'Add Stoichiometry' option. You can upload a file with a determined stoichiometry to be applied to the macrocomplex. As said before, the format of this file has to be the ID of the sequence, concordant with the one in the FASTA file, followed by ':' and a number. This number will be the number of times the corresponging sequence will be in the final complex.

<!-- ADD A GIF OF THE GUI OPERATION -->

Finally, in the top menu you can consult the ComplexBuilder 'Help' as well.

## Examples


## Limitations


## Team
| <a href="https://github.com/Paulagomis" target="_blank">**Paula Gomis Rosa**</a> | <a href="https://github.com/Argonvi" target="_blank">**Arturo González Vilanova**</a> | <a href="https://github.com/MartaLoBalastegui" target="_blank">**Marta López Balastegui**</a> |
| :---: |:---:| :---:|
| [![PaulaGomis](https://avatars2.githubusercontent.com/u/60719236?s=400&v=4)](https://github.com/Paulagomis)    | [![ArturoGonzalez](https://avatars1.githubusercontent.com/u/59646158?s=400&v=4)](https://github.com/Argonvi) | [![MartaLopez](https://avatars3.githubusercontent.com/u/44771228?s=400&v=4)](https://github.com/MartaLoBalastegui)  |

