<img src="https://cdn.pixabay.com/photo/2017/10/25/06/13/protein-icon-2887050_960_720.png" title="ComplexConstructorLogo" alt="ComplexConstructorLogo" height="100" width="100">

<!-- [![FVCproductions](https://avatars1.githubusercontent.com/u/4284691?v=3&s=200)](http://fvcproductions.com) -->

<!--***INSERT GRAPHIC HERE (include hyperlink in image)***-->

# Complex Constructor

> Generate macrocomplexes superimposing paired interacting elements

> protein - protein or protein - DNA

> the result is stored in a PDB file

## Table of Content

- [Why building complexes?](#whybuildingcomplexes?)
- [Installation](#installation)
- [Options](#options)
- [Examples](#examples)
- [Performance](#performance)
- [Limitations](#limitations)
- [Team](#team)


## Why building complexes?

The total number of proteins in humans is around 20K. From this quantity we already know the structure of 4K proteins. However, for 6K we have good templates and for other 6K of them, we have reasonably good templates; that means we have more or less 50-70% of the structure of human proteins covered.

On the other hand, the number of interactions that may occur between all the proteins is still unknown. There are some studies estimating 130-650K interactions, while others estimate more than millions. In this context, 120K interactions are confirmed from experiments and only 7K structures of complexes are known, which is a very small proportion of the total number of possible interactions anyway.

So, the complex structure coverage, including homologs, is around the 30%, much less than the mentioned coverage of monomers. That means there is still a long way to go into the complex-structure study field. 

The 3D structure of interacting proteins is necessary for them to develop their function, it is also fundamental for molecular recognition and contributes to the complexity of protein interaction networks too. However, identifying the structure of interacting proteins, complexes, is not an easy task. In monomers, homology modelling can be used to nearly cover all folds to determine a protein structure but, for interactions it is not enough. 

Complex Constructor tries to generate macrocomplex structures. To do so, the superimposing technique is used: it receives a list of pdb files, each of these files contains the structure of an interacting pair and, by superimposing the common elements of different pairs, it builds the final structure. Although there are other methods to generate macrocomplexes, the superimposition is fast and effective. Furthermore, not only protein-protein interacting pairs can be analyzed, but also proteins with DNA, to end up generating a macrocomplex structure of proteins and DNA chains. 

The core of Complex Constructor is the `constructor` function, which is the responsible of the building process. The pdbs of the interacting pairs are analyzed and classified using the information of the FASTA file. The common chains in different pdb files are identified: they will have more than 99% of identity in their sequences and an RMSD of less than 0.05 after superimposing their structures. After this classification, `constructor` begins to append elements to the macrocomplex: it starts with a pair, then it looks for another pair with a common element with the starting one and, if after superimposing the common element, the structure has not clashed, it appends the new element. Like this, the macrocomplex have now three elements. The next step will repeat the process looking for another pair, superimposing the common element and checking the possible clashes. The program is described [here](#performance).

In addition, as a lot of complexes follow a determined stoichiometry, Complex Constructor can also make the construction following it. If a stoichiometry is defined, the final structure will contain the exact number of elements indicated by it.

Finally, a new pdb file with the resulting structure of the macrocomplex is stored in the output folder called `outputFolder_model.pdb`.

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

### Install from pip

- ComplexConstructor can also be installed with `pip`

```shell
$ pip install complexconstructor
```
<!-- Change name of the repo! ComplexConstructor, pip ComplexConstructor?? -->
## Options

ComplexConstructor can be run using command-line arguments or using the graphical interface.

### Command-line

You can introduce the different arguments via command-line:

#### Mandatory arguments

To execute Complex Constructor three arguments are required:

- `-fa` `--fasta`: FASTA file with the sequences of the proteins or DNA that will conform the complex.

- `-pdb` `--pdbDir`: directory containing the pdb files with the structure of the pairs that will conform the complex.

- `-o` `--output`: directory name where the complex results will be stored. 

> Note that, if the output directory already exists, the results will be overwritten.


#### Optional arguments

- `-v` `--verbose`: show the detailed progression of the building process in a file called 'ComplexConstructor.log'. It will be stored in the output directory.

- `-st` `--stoichiometry`: File in .txt format containing the stoichiometry of the complex. The information of the stoichiometry must be in the following format: the ID of the sequence chain (concordant with the FASTA file ID) followed by ' : ' and the number of times it has to be present in the complex.

     `ID_FASTA_file` : `stoichiometry` One per line . Take a look at some examples [here](#examples).
     In case that this option is not used, the program will build the complex using each chain only once. 
### Graphical interface

Otherwise, the macrocomplex can also be built using the graphical interface:

```shell
$ python3 complexconstructor -gui
```

>In this case just the `-gui` tag is needed!

- To build the macrocomplex fill in the main window requirements. As for running the program via command-line, a FASTA file with the sequences and a directory with the pdb files of the interacting pairs are required. In addition, a name for the folder where the results will be stored is needed, after typing it you should confirm it. 

> Note that, to select the PDB directory you have to enter in the desired directory and then select it. You can see how it works in the [examples](#examples).

Furthermore, additional options can be set:

<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/gui.gif" alt="GUI" style="max-width:70%;" >

- In the main window you can specify if you want to create a log file where the process of the execution will be displayed. It will be stored in the output directory.

- In the top menu, in the 'Options' dropdown, there is the 'Add Stoichiometry' option. You can upload a file with a determined stoichiometry to be applied to the macrocomplex. As said before, the format of this file has to be the ID of the sequence, concordant with the one in the FASTA file, followed by ':' and a number. This number will be the number of times the corresponging sequence will be in the final complex.

Finally, in the top menu you can consult the Complex Constructor 'Help' as well.

## Examples

As said before, to generate any macrocomplex structure it is required the FASTA file and a direcory with paired structures in pdb format. All the data needed to execute the following examples is already located in the folder `examples`. Tu run them, we need to be in the folder of the package. 

```shell
$ cd ComplexConstructor
```

You can make sure you are in the correct folder with `ls` command.

```shell
$ ls
```

If you obtain the following directories: `assets`,  `complexconstructor`,  `examples`, together with other items, you are in the correct place!

> Inside the folder `examples` there are all the examples we will describe in this tutorial in the corresponding directories: `1gzx`, `5fj8`, etc. Each of them contains the required files to run Complex Constructor and also a pdb file of the reference structure, `exampleRef.pdb`, in order to check the superimposition between the constructed model and the real structure. 
<!--If Chimera or ICM is used, maybe the Ref files should be in the output foler??? -->

### 1GZX

Let's begin with the first example, the protein 1GZX. It is a small complex composed by two different aminoacid chains and each one of them appears two times (stoichiometry 2A2B), so the final structure has four chains. In order to perform the construction of T-state haemoglobin (1GZX), we need a .txt file where the stoichiometry is explicited. This file can be found in the folder `examples/1gzx` and it is called `1gzx_st.txt`:

```shell
1GZXa:2
1GZXb:2
```

1. Command line execution:

```shell
$ python3 complexconstructor -fa examples/1gzx/1gzx.fa -pdb examples/1gzx/1gzxDir -o 1GZX -st examples/1gzx/1gzx_st.txt -v
```

- `-fa`, mandatory: followed by the location of the FASTA file `1gzx.fa`.

> This file contain two IDs followed by the corresponding sequence, e.g. `1GZXa`, `1GZXb`. Note that, the IDs in `1gzx_st.txt` have to be concordant with them.

- `-pdb`, mandatory: followed by the location of directory with paired structures in pdb `1gzxDir`.

> In this case inside this folder we should have at least three pdb files, e.g. `1gzx_AB.pdb`, `1gzx_AC.pdb`, `1gzx_AD.pdb`. If there are redundant pairs they won't be considered. 

- `-o`, mandatory: followed by the name given to the output directory where the results will be stored, `1GZX`. 

> If a directory with the same name already exists in the working folder it will be overwritten.

- `-st`: followed by the stoychiometry information of the complex 2A2B that is in the file `1gzx_st.txt`.

- `-v`: turn ON the the verbose option. It is always recommended to create a logfile where the process information will be displayed. To deactivate the creation of the logfile, do not add the `-v` flag. 

2. Graphical interface execution:

<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/1gzxExample/1gzx.gif" title="1gzxGUI" alt="1gzxGUI" style="max-width:60%;" >

The resulting structure is stored in the directory `1GZX`, file `1GZX_model.pdb`.

| **Complex Constructor** | **Reference structure** | **Superimposition** |
| :---: | :---: | :---: |
|<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/1gzxExample/1gzxCC.png" title="1gzxCC" alt="1gzxCC" >|<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/1gzxExample/1gzxREF.png" title="1gzxREF" alt="1gzxREF" >|<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/1gzxExample/1gzxREF_CC.png" title="1gzxREF_CC" alt="1gzxREF_CC" style="max-width:92%;">


We can observe that the resulting structure from Complex Constructor fits the reference downloaded from PDB quite well. The RMSD of the second chains of both model and reference, computed with ICM after supeimposing the first chains, is zero. 

### 3KUY

3KUY is a complex composed by a DNA coil and a core made of protein chains. There are four aminoacid chains and one of nuclotides, all with stoichiometry two, which make a total of 10 chains. The procedure to run this example is the same as before. The data to construct the complex is inside the folder `examples/3kuy`. Execution with command-line arguments:

```shell
$ python3 complexconstructor -fa examples/3kuy/3kuy.fa -pdb examples/3kuy/3kuyDir -o 3KUY -st examples/3kuy/3kuy_st.txt -v
```

The resulting structure in directory `3KUY`, file `3KUY_model.pdb`.

> To execute this example, and the same for all the rest, with the graphical interface repeat the same process as for the previous example with the inputs of the folder `examples/3kuy`.

| **Complex Constructor** | **Reference structure** | **Superimposition** |
| :---: | :---: | :---: |
|<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/3kuyExample/3kuyCC.png" title="3kuyCC" alt="3kuyCC" >|<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/3kuyExample/3kuyRef.png" title="3kuyref" alt="3kuyRef" >|<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/3kuyExample/3kuyRef_CC.png" title="3kuyRef_CC" alt="3kuyRef_CC" style="max-width:92%;">

The whole complex is correctly constructed. 
<!-- add info here, Zscore, energy...-->

### 4R3O

A bigger case is the complex 4R3O, but just made of aminoacid chains. It is also symmetric and composed by 14 chains, stoichimetry two for all of them, a total of 28 chains in the complex. Its input data in `examples/4r3o`.

```shell
$ python3 complexconstructor -fa examples/4r3o/4r3o.fa -pdb examples/4r3o/4r3oDir -o 4R3O -st examples/4r3o/4r3o_st.txt -v
```

The resulting structure in directory `4R3O`, file `4R3O_model.pdb`.

| **Complex Constructor** | **Reference structure** | **Superimposition** |
| :---: | :---: | :---: |
|<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/4r3oExample/4r3oCC.png" title="4r3oCC" alt="4r3oCC" >|<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/4r3oExample/4r3oRef.png" title="4r3oref" alt="4r3oRef" >|<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/4r3oExample/4r3oRef_CC.png" title="4r3oRef_CC" alt="4r3oRef_CC" style="max-width:92%;">

In this case, as well, the constructed model fits the reference structure. <!-- add info here, Zscore, energy...-->

### 5FJ8

The next complex is composed by aminoacid and nucleotide sequences but, in this case, it is not symmetric. It is composed by 20 chains, all just present once in the structure. The required inputs to the construction are in `examples/5fj8`. 

```shell
$ python3 complexconstructor -fa examples/5fj8/5fj8.fa -pdb examples/5fj8/5fj8Dir -o 5FJ8 -st examples/5fj8/5fj8_st.txt -v
```

The resulting structure in directory `5FJ8`, file `5FJ8_model.pdb`.

| **Complex Constructor** | **Reference structure** | **Superimposition** |
| :---: | :---: | :---: |
|<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/5fj8Example/5fj8CC.png" title="5fj8CC" alt="5fj8CC" >|<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/5fj8Example/5fj8Ref.png" alt="5fj8Ref" >|<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/5fj8Example/5fj8Ref_CC.png" title="5fj8Ref_CC" alt="5fj8Ref_CC" style="max-width:92%;">

The model and the reference are superimposed with very little differences between both structures. In this particular case, the aminoacid chain Q, had several aminoacids labelled as 'unknown' in the pdb files and as 'X' in the FASTA sequence. To deal with this, we had to take out this aminoacids from input files, and so, the Q chain is partly constructed in the model. Nevertheless, the rest of the structure is correctly reproduced.
<!-- add info here, Zscore, energy...-->

### 6GMH
Another non-symetric example, but also with DNA sequences, is 6GMH. It has 17 aminoacid chains and 3 DNA chains.

```shell
$ python3 complexconstructor -fa examples/6gmh/6gmh.fa -pdb examples/6gmh/6gmhDir -o 6GMH -st examples/6gmh/6gmh_st.txt -v
```

Resulting structure in directory `6GMH`, file `6GMH_model.pdb`.

| **Complex Constructor** | **Reference structure** | **Superimposition** |
| :---: | :---: | :---: |
|<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/6gmhExample/6gmhCC.png" title="6gmhCC" alt="6gmhCC" >|<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/6gmhExample/6gmhREF.png" title="6gmhref" alt="6gmhRef" >|<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/6gmhExample/6gmhREF_CC.png" title="6gmhRef_CC" alt="6gmhRef_CC" style="max-width:92%;">

In this case, the pdb strucure had bigger regions labelled as 'unknown'. As those regions are not recognized by Complex Constructor, we removed them from the input pdb pairs and from the FASTA sequences, as in the previous example. That is why the 'Q' chain is not completelly constructed: from the original number of 884 aminoacids, after the deletion of the 'unknowns', 584 remained. But, here as well, the correctly characterized aminoacids, and all the DNA sequences, are constructed in the right way as we can see in the image of the superimposition. 
<!-- add info here, Zscore, energy...-->

## Performance

<img src="https://raw.githubusercontent.com/Argonvi/SBI-PYT_Project/master/assets/ComplexBuilderDiagram.jpg" title="ComplexConstructorLogo" alt="ComplexConstructorDiagram" >

- Before adding a new chain to the macrocomplex the number of clashes between the new chain and the previous structure is checked. The function `sequence_clashing` finds how many CA atoms from the new chain are closer than 2 angstroms to any other CA atom of the previous macrocomplex, this is, the number of clashes. If the number of clashes is above 20, the new chain won't be added to the macrocomplex. 

### Structure of the package

#### ComplexConstructor

- `main.py`: main module of Comple Constructor. It connects with all the rest of modules. 

- `utilities.py`: module with all the functions needed to run the construction of macrocomeplexes, as `constructor` or `superimpositor`.

- `argparser.py`: reads and organizes the command-line arguments.

- `interface.py`: contains the configuration of the graphical interface and passes the inputs to the main module.

- `logProgress.py`: this module creates the log file and its contents when the verbose option is selected.

- `helpText.py`: it is a module with the help text to be added to the graphical interface 'Help' option.


## Limitations

- The main limitation of the program is the preparation of the input files. The FASTA file must have each ID just once, the repited sequences have to be erased. Furthermore, if the pairs of the pdb files are not correctly separated the program won't be able to identify the IDs of the FASTA file in the pdb files. Finally, for each ID, the sequences of the FASTA file must be highly similar to the one of the pdb, otherwise they won't be identified as the same.

- The aminoacids labelled as 'unknown' in the pdb files and as 'X' in the fasta sequences, cannot be correctly identified so they should not be present in the input files. As a consequence, those chains are just partly constructed in the final model.

- Although the similarity of chains from different pdb files is double checked, only the ones with sequence identity over 99%, and also with an RMSD under 0.05 when both structures are superimposed, the small differences between them may be propagated when building big structures. It could be partly solved if the similarity conditions were harder or with small structure corrections everytime a chain is added.

- The refinement of the structures obtained would have been also very usefull. 

- We would have also liked to check the Z-scores of the resulting structures through the energy analisys to check the validity of the results.


## Team
| <a href="https://github.com/Paulagomis" target="_blank">**Paula Gomis Rosa**</a> | <a href="https://github.com/Argonvi" target="_blank">**Arturo González Vilanova**</a> | <a href="https://github.com/MartaLoBalastegui" target="_blank">**Marta López Balastegui**</a> |
| :---: |:---:| :---:|
| [![PaulaGomis](https://avatars2.githubusercontent.com/u/60719236?s=400&v=4)](https://github.com/Paulagomis)    | [![ArturoGonzalez](https://avatars1.githubusercontent.com/u/59646158?s=400&v=4)](https://github.com/Argonvi) | [![MartaLopez](https://avatars3.githubusercontent.com/u/44771228?s=400&v=4)](https://github.com/MartaLoBalastegui)  |


