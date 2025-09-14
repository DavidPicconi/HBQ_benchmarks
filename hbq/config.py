from .parameters import load_parameters, load_geometry_data

# Compute the parameters of the diabatic PESs
V1D, grad, Hess, c3, c4, W12, W13, W23 = load_parameters()

# Get geometrical data
n_atoms, geo_S0, mass_tot, masses, TR, C, Cm1, Omega, Ueff, w0 = load_geometry_data()

# Folder for the parameter data
import os.path
parameters_dir = os.path.join(os.path.dirname(__file__), 'surface_parameters')