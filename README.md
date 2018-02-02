[![Build Status](https://travis-ci.com/CINPLA/KNPsim.svg?token=CMEYLVk7cYLKvFcNYS19&branch=master)](https://travis-ci.com/CINPLA/KNPsim)

## KNPsim
KNPsim is a implementation of the KNP framework (see Halnes et al., https://www.ncbi.nlm.nih.gov/pubmed/27820827), using the FEniCS software for finite element problems.
KNP is a framework for solving ion concentration dynamics by assuming the Nernst-Planck equation and electroneutrality. Neuronal input is provided as point sources to the system.

# Installation
Work in progress:
```
conda create -n knpsim -c conda-forge fenics h5py matplotlib python=3
source activate knpsim
python setup.py install
```

# Quick introduction
TODO

# Examples
Examples of usage is provided in the folder `src/simulations`.

# References
TODO
