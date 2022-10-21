In this folder results of the simulations, related parameters and notebooks to perform the analyses are reported.

Each subfolder corresponds to a given combination of parameters for which MD simulations were executed. In particular, simulations consisted in relaxation of the initial structure already minimized with AMBER. Each folder contains files bonds.pkl, angles.pkl, dihedrals.pkl and rmsd.npy, that store respectively the length of bonds and the amplitude of angles and dihedrals for each combination of particles, and the RMSD with respect to the initial RNA molecule. All these files are generated with the function main() in relaxation_analysis.ipynb, after giving the correct folder as input. Some folders also contain a subfolder Relaxed_conformations, with the final structures obtained after the relaxation.

The methods used are the following:

- [Amber](Amber): original minimized structures;
- [HiRE_old](HiRE_old): initial values of parameters;
- [FreeMin](FreeMin): free minimization, with learning rate equal for all parameters;
- [FreeMin_newangles](FreeMin_newangles): free minimization, but with large initial values for angles involving purines;
- [PSGD](PSGD): preconditioned stochastic gradient descent;
- [NoGeom](NoGeom): fixed geometric parameters;
- [Less_dof](Less_dof): global parameters fixed to 1 to reduce degeneracy;
- [Tanh](Tanh): little variation allowed for geometric parameters, following a tanh function.

The other notebooks are used to plot results for the RMSD, analyse AMBER distributions in terms of Gaussian fitting and compare results for the different FreeMin implementations.

Finally, [MD_parameters](MD_parameters) contains files print_table.py and parameters_from_txt.py to generate parameters tables from .pth files and substitute them in the inputs for the HiRE model, and [MD_utils](MD_utils) contains subfolders with parameters_table.txt used for each simulation.  
