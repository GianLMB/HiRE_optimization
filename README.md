## Introduction

This project aims at optimizing the parameters of a coarse-grained force field with methods based on machine learning (ML) and physical intuition. The HiRE-RNA force field [1], developed at Université Paris Cité (add github repository), is designed to reproduce the experimental results of RNA folding. It contains specific functional forms that describe anergy and forces both local and non-bonded interactions, with particular attention for stacking and hydrogen bonds in order to take into account also non-canonical base pairings. Overall, involved parameters are more than 300, related to energy couplings, equilibrium values and particle-specific coefficients.
In this first phase of the project I focussed on parameters related to local interactions, the optimization of which is based on the matching between HiRE and AMBER [2] energies. I designed the machine learning setup from scratch, including:

1. Creation and storage of a suitable dataset to train the model, with the extraction of relevant features from PDB files of RNA sequences (coordinates and indexes of particles, type of interactions, energy labels, etc.);
2. Translation of the HiRE-RNA code from Fortran to Python, to exploit PyTorch's automatic differentiation engine;
3. Definition of the ML models to compute the energy with HiRE for parameters optimization;
4. Analysis of results through molecular dynamics (MD) simulations.

The language used is Python 3.8 and the ML model is written using Pytorch library. For the dataset creation and analysis of results, [biopython](https://biopython.org/) and [MDtraj](https://www.mdtraj.org/1.9.8.dev0/index.html) packages are required.

## Directories

[Models](Models) contains the codes and data for different implementations of ML models, based on the optimized parameters and the techniques used.

[MD_simulations](MD_simulations) contains the results and analyses for MD simulations, with different choice of parameters.

[Structure_DB_Amber_HiRE](Structure_DB_Amber_HiRE) contains the PDB files of RNA structures used to construct the dataset. In particular it is divided into 2 subfolders:

1. [FA_PDB](Structure_DB_Amber_HiRE/FA_PDB), with original RNA sequences, the subsequences obtained parsing PDB files (Cropped_PDB) and the corresponding energies computed with AMBER (Cropped_energies);
2. [Prep_pureHire](Structure_DB_Amber_HiRE/Prep_pureHire), with the script to generate coarse-grained structures for each sequence and eventually compute energies with HiRE, that are then stored in homonymous directories in Cropped_Output.

[HiRE_lib](HiRE_lib) contains the source code in Fortran and the instructions to build the HiRE model.


## References

[[1]](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00200)
Tristan Cragnolini, Yoann Laurin, Philippe Derreumaux, and Samuela Pasquali.
Coarse-Grained HiRE-RNA Model for ab Initio RNA Folding beyond Simple Molecules, Including Noncanonical and Multiple Base Pairings.
Journal of Chemical Theory and Computation,
11(7):3510–3522, 2015.

[[2]](https://ambermd.org/index.php) David Case, Ido Ben-Shalom, S.R. Brozell, D.S. Cerutti, Thomas Cheatham,Thomas Darden, Robert Duke, Delaram Ghoreishi, Michael Gilson, H. Gohlke, Andreas Gotz, D. Greene, Robert Harris, N. Homeyer, Yandong Huang, Saeed Izadi, Andriy Kovalenko, Tom Kurtzman, and P.A. Kollman.
Amber 2022.
Technical report, University of California, 2022.
