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

- Background theory


## Installation

- All the `code` required to get started
- Images of what it should look like

### Clone

- Clone this repo to your local machine using `https://github.com/fvcproductions/SOMEREPO`

### Setup

- If you want more syntax highlighting, format your code like this:

> update and install BioPython 

```shell
$ brew update
$ brew install fvcproductions
```



## Options

- ComplexBuilder can be run using command-line argument or using the graphical interface.

### Graphical interface

- To build the macrocomplex with the graphical interface 
```shell
$ python3 ComplexBuilder.py -gui
```
- In this case just the -gui tag is needed!

- To build the macrocomplex fill in the main window requirements. As for running via command-line, a FASTA file with the sequences and a directory with the PDB files of the interacting are required. In addition, a name for the folder where the results will be stored is needed, after typing it you should confirm it. 

- Note that, as seeing in the following demonstration, to select the PDB directory you have to enter in the desired directory and then select it. 

- Furthermore, additional options can be set:

> In the main window you can specify if you want to create a log file where the process of the execution will be displayed. 

> In the top menu, in the 'Options' dropdown there is the 'Add Stoichiometry' option. You can upload a file with a determined stoichiometry to be applied to the macrocomplex. The format of this file has to be the ID of the sequence, concordant with the one in the FASTA file, followed by ':' and a number. This number will be the number of times the corresponging sequence will be in the final complex. See an example [here](#optional arguments).

<!-- ADD A GIF OF THE GUI OPERATION -->

- Finally, in the top menu you can consult the ComplexBuilder 'Help' as well.

### Command-line

- Otherwise you can introduce the arguments via command line.

#### Mandatory arguments

- To execute ComplexBuilder three arguments are required:

> -fa: FASTA file with the sequences of the proteins or DNA that will conform the complex.

> -pdb: diretory containing the PDB files with the structure of the pairs that will conform the complex.

> -o: directory name where the complex results will be stored. 

#### Optional arguments

> -v: show the detailed progression of the building process in a file called 'ComplexBuilder.log'.

> -st: File containing a determined stoichiometry to the complex. The information of the stoichiometry must be: the ID of the sequence chain (concordant with the FASTA file ID) followed by the number of times it has to be present in the complex after ':'
ID_as_FASTA_file : stoichiometry (one per line) in format .txt. Example for a stoichiometry 2A2B for '1GZX':

```shell
1GZXA:2
1GZXB:2
```

## Examples


## Limitations


## Team
| <a href="https://github.com/Paulagomis" target="_blank">**Paula Gomis Rosa**</a> | <a href="https://github.com/Argonvi" target="_blank">**Arturo González Vilanova**</a> | <a href="https://github.com/MartaLoBalastegui" target="_blank">**Marta López Balastegui**</a> |
| :---: |:---:| :---:|
| [![PaulaGomis](https://avatars2.githubusercontent.com/u/60719236?s=400&v=4)](https://github.com/Paulagomis)    | [![ArturoGonzalez](https://avatars1.githubusercontent.com/u/59646158?s=400&v=4)](https://github.com/Argonvi) | [![MartaLopez](https://avatars3.githubusercontent.com/u/44771228?s=400&v=4)](https://github.com/MartaLoBalastegui)  |

