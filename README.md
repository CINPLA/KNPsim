[![Build Status](https://travis-ci.com/CINPLA/KNPsim.svg?token=CMEYLVk7cYLKvFcNYS19&branch=master)](https://travis-ci.com/CINPLA/KNPsim)

## KNPsim
KNPsim is a implementation of the KNP framework (see Halnes et al., https://www.ncbi.nlm.nih.gov/pubmed/27820827), using the FEniCS software for finite element problems.
KNP is a framework for solving ion concentration dynamics by assuming the Nernst-Planck equation and electroneutrality. Neuronal input is provided as point sources to the system. This code was used to produce the results
presented in (TODO)

# Installation
Work in progress:
```
conda create -n knpsim -c conda-forge fenics ipython h5py matplotlib python=3
source activate knpsim
python setup.py install
```

# Testing
Run the command `pytest` from the top directory. Make sure you have `pytest`
installed (`conda install pytest`).

# Reproducing the results of the paper
In order to reproduce the results of (TODO: future paper), run the script `run_all_examples_and_gather_data.sh`. Then go to `KNPsim/make_the_figures`,
open the notebooks and run all the cells in each.
Note that All the figures have been formatted in Inkscape for publication
purposes.

# How to set up an electrodiffusion problem in knpsim
A script using knpsim will usually go through the following steps:
1. Create an instance of `Geometry` (will involve making or loading a mesh).
2. Create an instance of `Simulator`.
3. Create instances of the ion species.
4. (optional) Set up current point sources
5. Create an instance of `Time_solver`
6. Initialize the simulation by calling `simulator.initialize_simulator()`
7. (optional) Create an instance of `State_saver`, to save the result.
8. Run the simulation by calling `time_solver.solve()`
