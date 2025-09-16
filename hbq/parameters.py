"""
Module to read and process the parameter files.

@author: David Picconi
"""
import os.path
import numpy as np

###############################
def load_geometry_data():
    mass = {'H': 1.00782503223,
            'C': 12.0,
            'N': 14.00307400443,
            'O': 15.99491461957}               # https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
    
    DALTON_TO_AU = (1.0 / 5.485799090441e-4)   # https://physics.nist.gov/cgi-bin/cuu/Value?meu|search_for=electron+mass
    ANG_TO_BOHR  = (1.0 / 0.529177210544)      # https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0|search_for=bohr
    
    datapath = os.path.dirname(__file__)
    
    with open(os.path.join(datapath, 'geo_S0min.dat'), 'r') as fh:
        geom   = fh.readlines()
        geo_S0 = np.array([[float(x) for x in l.split()[1:]] for l in geom])
        masses = np.array([mass[l.split()[0]] for l in geom])
            
    n_atoms  = len(geo_S0)
    mass_tot = np.sum(masses)
    
    # Align to the center of mass
    rcm   = (masses @ geo_S0) / mass_tot
    geo_S0 = geo_S0 - rcm
    geo_S0 = geo_S0.reshape(n_atoms * 3) * ANG_TO_BOHR
    
    # Read the normal modes, their frequencies and the effective modes
    U0   = np.loadtxt(os.path.join(datapath, 'U0.dat'))   # normal modes
    Ueff = np.loadtxt(os.path.join(datapath, 'Ueff.dat')) # effective modes
    w0   = np.loadtxt(os.path.join(datapath, 'w0.dat'))   # S0 frequencies
    TR   = np.loadtxt(os.path.join(datapath, 'TR.dat'))   # translations/rotations
    
    Omega = Ueff.T @ np.diag(w0) @ Ueff
    
    masses_3 = np.c_[masses,masses,masses].reshape(n_atoms * 3) * DALTON_TO_AU
    sqrt_masses = np.sqrt(masses_3) 
    TR = np.diag(sqrt_masses) @ TR[:,[3,4,5]]
    
    C   = Ueff.T @ np.diag(np.sqrt(w0)) @ U0.T @ np.diag(sqrt_masses)
    Cm1 = np.diag(1 / sqrt_masses) @ U0 @ np.diag(1 / np.sqrt(w0)) @ Ueff
    
    return n_atoms, geo_S0, \
           mass_tot * DALTON_TO_AU, masses * DALTON_TO_AU, \
           TR, C, Cm1, Omega, Ueff, w0


###############################
def interpolate(scanFile):
    """
    Read the data from datafile and perform a quadratic interpolation
    """    
    datapath = os.path.join(os.path.dirname(__file__), 'surface_parameters')
    
    from numpy import vectorize
    if not os.path.isfile(os.path.join(datapath, scanFile)):
        return vectorize(lambda x: 0.0)
    
    with open(os.path.join(datapath, scanFile), 'r') as FileTmp:
        scan = FileTmp.readlines()
    
    Q1 = [float(x.split()[0]) for x in scan]        
    V  = [float(x.split()[1]) for x in scan]
    
    from scipy.interpolate import interp1d
    
    return interp1d(Q1, V, kind = 'quadratic',
                    fill_value = 'extrapolate')
    
    
###############################
def load_parameters():
    """
    Read and interpolate all the parameters
    """
    nAtoms = 24
    nModes = 3 * nAtoms - 6
    
    #####
    # One-dimensional cuts
    V1D = [0] * 4
    for iS in range(4):
        scanFile = 'V' + str(iS) + '.dat'
        V1D[iS] = interpolate(scanFile)
        
    #####
    # Gradients
    grad = [[0 for iM in range(nModes)] for iS in range(4)]
    for iS in range(4):
        for iM in range(nModes - 1):
            scanFile = 'Grad_' + str(iS) + '_' + str(iM + 1) + '.dat'
            grad[iS][iM + 1] = interpolate(scanFile)
    
    #####
    # Hessians
    Hess = [[[0 for iM in range(nModes)] for jM in range(nModes)]\
              for iS in range(4)]  
    for iS in range(4):
        for iM in range(nModes - 1):
            for jM in range(iM, nModes - 1):
                scanFile = 'Hess_' + str(iS) \
                           + '_' + str(iM + 1) + '_' + str(jM + 1) + '.dat'
                Hess[iS][iM + 1][jM + 1] = interpolate(scanFile)
                
    #####
    # Cubic and quartic terms
    c3 = [0 for iS in range(4)]
    c4 = [0 for iS in range(4)]
    for iS in range(4):
        scanFile = 'c3_' + str(iS) + '.dat'
        c3[iS] = interpolate(scanFile)
        scanFile = 'c4_' + str(iS) + '.dat'
        c4[iS] = interpolate(scanFile)
        
    #####
    # Diabatic couplings
    W12 = interpolate('coup12.dat')
    W13 = [interpolate('coup13_1.dat'),
           interpolate('coup13_2.dat'),
           interpolate('coup13_3.dat'),
           interpolate('coup13_4.dat')]
    W23 = [interpolate('coup23_1.dat'),
           interpolate('coup23_2.dat'),
           interpolate('coup23_3.dat'),
           interpolate('coup23_4.dat')]
    
    ############################
    return V1D, grad, Hess, c3, c4, W12, W13, W23
    