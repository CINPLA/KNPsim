#!/usr/bin/env bash

# KNP long
python step_concentration.py --relax 0.9

# PNP long
python step_concentration.py -p poisson --start 0 --end 100e-6

# Nofield long
python step_concentration.py -p zero --start 0 --end 100e-6 -m False --relax 1 --rtol 1e-3

# KNP zoom
python step_concentration.py -z True --dt 1e-10 -T 1e-7 --start 49.9e-6 --end 50.1e-6 --rtol 1e-5

# PNP zoom
python step_concentration.py -z True -p poisson --dt 1e-10 -T 1e-7 --start 49.9e-6 --end 50.1e-6 --rtol 1e-5

# Modified diffusion
python step_concentration.py -p zero --start 0 --end 100e-6 --relax 1 -m True --rtol 1e-6
