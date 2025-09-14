# HBQ Benchmarks

This repository provides Python routines to compute the diabatic electronic Hamiltonian of **HBQ** (2-(2'-hydroxyphenyl)benzothiazole), designed for benchmarking nonadiabatic dynamics methods.  
The model was constructed as described in *Photochem. Photobiol. Sci.* **20**, 1455–1473 (2021).  
[https://doi.org/10.1007/s43630-021-00112-z](https://doi.org/10.1007/s43630-021-00112-z)

## Setup
To set up the package:
1. Download the repository.  
2. Extract the surface parameter files:  
```
tar -xvf surface_parameters.tar
```

3. The final folder structure should look like:
```
hbq/
├── surface_parameters/
│   ├── c3_0.dat
│   ├── c3_1.dat
│   └── ...
├── __init__.py
├── config.py
├── hamiltonian.py
└── ...
```

## Usage
You can import the hbq package directly into your Python code.

For step-by-step examples of setting up calculations and running the routines, see the Jupyter notebooks in the Examples/  folder.

## Citation
If you use this package in your research, please cite:

D. Picconi, Photochem. Photobiol. Sci. 20, 1455–1473 (2021).
https://doi.org/10.1007/s43630-021-00112-z
