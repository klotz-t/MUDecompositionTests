### Description
This repository contains Matlab codes to replicate the data presented in *Klotz and Rohlen (2025) "Revisiting convolutive blind source separation for identifying spiking motor neuron activity: From theory to practice".* 


### Installation 
1. You need a working MATLAB installation

2. Clone the repository:
```bash
git clone https://github.com/klotz-t/MUDecompositionTests.git
```

### Requirements
1. Matlab2024a
2. Signal Processing Toolbox
3. Statistics and Machine Learning Toolbox
4. Parallel Computing Toolbox

### Repository structure
```
MUDecompositionTests/
├── experimental-simulation/     % Contains the scripts used to extract the experimental MUAPs (+ a copy of the MUAPs we have used)
├── Figures/                     % Scripts for replicating the simulations and figures presented in the paper (+ replication data)
├── Functions/                   % Functions used for the upper bound decompositions + utility functions
├── LIF Model/                   % Implementation of the utilized leaky-integrate and fire model and functions for simulating EMG signals
└── README.md
```

### Run an in silico experiment and decompose the data
1. Open Matlab and browse to the folder hosting a copy of the repository

2. Browse to the Figures folder by entering the following command into your terminal:
```bash
cd Figures/
```
3. Run a script of interest via your Matlab terminal, e.g.,:
```bash
generate_figure4
```

### DOI
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14824963.svg)](https://doi.org/10.5281/zenodo.14824963)


