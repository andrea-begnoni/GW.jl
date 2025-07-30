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

To download the repository

```
git clone https://github.com/andrea-begnoni/GW.jl.git
```

You are now ready to go. We suggest you to start from `easy_GWJulia.ipynb`

## Structure of the repository

Directories:

- `src`: where all the source code is located
- `test`: contains some test functions
- `useful_files`: where the sensitivities of the interferometers and some files needed for the waveforms are kept
- `catalogs`: automatically generated during the first run, used to store generated catalogs
- `output`: automatically generated during the first run, used to store the outputs (e.g. Fisher matrix, forecast uncertainties)

### The `src` directory

The `src` folder contains the 4 modules used by the code:

- `catalog.jl` contains population models and functions to generate and load the catalogs
- `detector.jl` contains all relevant info (e.g., positions, noises, patter functions) to output the signal measured by each detector (for the waveform uses `waveforms`, see below)
- `utils.jl` contains multi-purpose functions and the functions to analyze the Fisher Matrix
- `waveforms` folder containing all available waveforms models, tested to be compatible with their LAL counterparts
