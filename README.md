# GW 
[![Build Status](https://github.com/andrea-begnoni/GW.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/andrea-begnoni/GW.jl/actions/workflows/CI.yml?query=branch%3Amain)

![alt text](logo.png)
### Perform accurate and fast parameter estimation of Gravitational Waves sources using the Fisher Information Matrix (alias GWJulia).

Main features:
- Generate a catalog of compact binary sources
- Perform a Fisher Matrix analysis with automatic differentiation
- Analyze the results and obtain the errors on the parameters of the binary 
- Very fast, it takes less than a second to compute the Fisher Matrix, independently of the waveform model and the detector(s) configuration

## Where to start?

After the installation there is an extensive tutorial called `easy_GWJulia.ipynb`, we leave to it the description of all the features available in the code.


## Installation

Before proceding with the installation verify that you have Julia installed on your machine and we suggest to use the most recent stable version. You can find instructions on how to install Julia here: https://julialang.org/downloads/.

The code is also a julia package so it is very easy to install, enter in julia
```
git clone https://github.com/andrea-begnoni/GW.jl.git
```

You are now ready to go. We suggest you to start from `easy_GWJulia.ipynb`

## Structure of the repository

Directories:
- `src` where all the source code is located
- `catalogs` where the code saves the catalogs generated
- `output` where the outputs (e.g. Fisher matrix, expected errors) are saved
- `useful_files` where the sensitivities of the interferometers and some files needed for the waveforms are kept 


### The `src` directory

The `src` forlder contains the 4 modules that make the code:
- `waveform.jl` contains the waveforms models available, from LAL
- `detector.jl` it calls the waveform from the previous module and contains functions to do all the computations till the fisher matrix
- `catalog.jl` generates and reads the catalogs
- `utils.jl` contains multi-purpose functions and the fucntions to analyze the Fisher Matrix


