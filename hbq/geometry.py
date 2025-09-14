"""
Module to read and manipulate the geometry

@author: David Picconi
"""

import numpy as np

from .config import geo_S0, mass_tot, masses, C, Cm1, TR, Ueff

###############################


def align(geo_start, check_alignment = False):
    """
    Rotate the molecule so to eliminate the projection along
    the infinitesimal translational coordinates
    """
    geo = geo_start.copy()
    # Check the alignment
    if check_alignment:
        print('')
        print('ALIGN: Initial components on the rotational coordinates:')
        print(TR.T @ geo.flatten())
    
    max_iter = 200
    iC = 0
    eps = 1e-8
    for i in range(max_iter):
        ixyz = (i % 3)    # Cyclic rotations about the three axes
        if   ixyz == 0:
            geo_swap = np.array([-geo[:,2],  geo[:,1],  geo[:,0]]).T
        elif ixyz == 1:
            geo_swap = np.array([ geo[:,0],  geo[:,2], -geo[:,1]]).T
        elif ixyz == 2:
            geo_swap = np.array([ geo[:,1], -geo[:,0],  geo[:,2]]).T
        
        b = np.dot(TR[:,ixyz], geo.flatten())
        c = np.dot(TR[:,ixyz], geo_swap.flatten()) 
        if np.abs(b) < eps:
            iC += 1
        else:
            iC = 0
        if iC == 3: break
        
        theta = np.arctan(- b / c)
        s_rot = np.sin(theta)
        c_rot = np.cos(theta)
        geo = c_rot * geo + s_rot * geo_swap
        
    if i == max_iter - 1:
        print('WARNING: very distorted geometry.')
        print('It could not be completely aligned.')
        
    # Check alignment
    if check_alignment:
        print('ALIGN: %i iterations' % (i + 1))
        print('ALIGN: Final components on the rotational coordinates:')
        print(TR.T @ geo.flatten())        
    #
    return geo.flatten()



def cart_to_q(geo, check_alignment = False):
    """
    Converts the geometry from Cartesian coordinates into dimensionless effective modes
    """
    # Check the shape
    if geo.shape != (24,3):
        print('ERROR: the input geometry should be a (24,3) numpy array')
        return None
    
    # Center the center of mass to the origin
    rcm = (masses @ geo) / mass_tot
    geo_shift = geo - rcm
           
    # Align the molecule
    geo_shift_aligned = align(geo_shift, check_alignment)
    
    # Evaluate the displacement
    Q = C @ (geo_shift_aligned - geo_S0)
    
    #
    return Q



def q_to_cart(Q):
    """
    Converts the geometry from dimensionless effective modes into Cartesian coordinates
    """
    # Check the length
    if Q.shape != (66,):
        print('ERROR: the input should be a (66,) numpy array')
        return None
    
    geo_cart = geo_S0 + Cm1 @ Q
    
    #
    return geo_cart.reshape((24,3))



def nm_to_q(Q_nm):
    return Ueff.T @ Q_nm


def q_to_nm(Q):
    return Ueff @ Q
        

def nm_to_cart(Q_nm):
    Q = nm_to_q(Q_nm)
    return q_to_cart(Q)


def cart_to_nm(geo, check_alignment = False):
    Q = cart_to_q(geo, check_alignment)
    return q_to_nm(Q)
    