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

# How to set up an electrodiffusion problem in knpsim
A script using knpsim will usually go through the following steps:
1. Create an instance of Geometry (will involve making or loading a mesh).
2. Create an instance of Simulator.
3. Create instances of the ion species.
4. 

# Reproducing the results of the paper
In order to reproduce the results of (TODO: future paper), perform the
following steps:
1. Go to each subfolder in `KNPsim/examples` and run the examples. Each example
  has a README.md with details on how to run them.
2. Move the resulting h5 files to `KNPsim/make_the_figures/data`.
3. In `KNPsim/make_the_figures` there is a notbook for each figure. Open each
  notebook and run all the cells.
Note that All the figures have been formatted in Inkscape for publication
porposes.
# References
TODO
