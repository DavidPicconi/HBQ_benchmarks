"""
Module to evaluate the diabatic potential energy surface
for the HBQ molecule. 
The ground state and the excited states, as well as the diabatic couplings
between the excited states can be computed.

@author: David Picconi
"""

#from .geometry import geoS0, masses, Omega   # Load the S0 minimum geometry
                                     # and the atomic masses

# Initialise the PES parameters
from .config import V1D, grad, Hess, c3, c4, W12, W13, W23, Omega, parameters_dir

##################################
# Functions to evaluate the diabatic potentials 

def electronic_hamiltonian_q(Q):
    """
    Calculate the diabatic electronic Hamiltonian
    using the dimensionless normal modes as input.
    
    Input:
        Q  : The vector of dimensionless coordinates.
             It should have 66 dimensions
    
    Output
        W: The ground and excited state 4x4 diabatic potential matrix
           in hartrees.
           The couplings between the ground and the excited states are zero
    """
    
    n_modes = 66
    # Check the size
    if len(Q) != n_modes:
        print('ERROR: an array of lenght 66 was expected')
        return None
    
    # Potentials
    import numpy as np
    
    W = np.zeros((4,4))
    Q1 = Q[0]
    
    # diagonal potentials
    for i in range(4):
        V = V1D[i](Q1)

        # linear and quadratic terms
        for r in range(n_modes - 1):
            V += grad[i][r + 1](Q1) * Q[r + 1] \
               + 0.5 * Hess[i][r + 1][r + 1](Q1) * Q[r + 1]**2
            for s in range(r + 1, n_modes - 1):
                V += Hess[i][r + 1][s + 1](Q1) * Q[r + 1] * Q[s + 1]
                
        # cubic and quartic terms
        V += c3[i](Q1) * Q[8]**3 + c4[i](Q1) * Q[8]**4
        
        #
        W[i,i] = float(V)
        
    # diabatic couplings
    W[1,2] = W12(Q1)
    
    W[1,3] = W13[0](Q1) * Q[1] \
           + W13[1](Q1) * Q[2] \
           + W13[2](Q1) * Q[3] \
           + W13[3](Q1) * Q[4]
           
    W[2,3] = W23[0](Q1) * Q[1] \
           + W23[1](Q1) * Q[2] \
           + W23[2](Q1) * Q[3] \
           + W23[3](Q1) * Q[4]
           
    W[2,1] = W[1,2]
    W[3,1] = W[1,3]
    W[3,2] = W[2,3]
    
    #
    return W


    
def electronic_hamiltonian_cart(geo, check_alignment = False):
    """
    Calculate the diabatic potentials using
    Cartesian coordinates (in bohr) as input.
    
    Input:
     geo : The geometry of HBQ in bohr.
             It should have shape (24,3)
     check_alignment : if True a message is printed to check that
                      the molecule was rotated correctly
    
    Output
        W: The ground and excited state 4x4 diabatic potential matrix
           in hartrees.
           The couplings between the ground and the excited states are zero
    """
    
    from .geometry import cart_to_q
    Q = cart_to_q(geo, check_alignment)
    #
    return electronic_hamiltonian_q(Q)



def electronic_hamiltonian_nm(Q_nm):
    """
    Calculate the diabatic potentials using
    dimensionless normal modes as input.

   Input:
    Q_nm : The normal modes displacement vector (dimensionless)
           It should be a 66-dimensional numpy array
    """
    
    from .geometry import nm_to_q
    Q = nm_to_q(Q_nm)
    #
    return electronic_hamiltonian_q(Q)



#################################################
# Function to generate an operator file for MCTDH

def mctdh_operator(modes, op_file = 'hbq.op', parameters_path = None):
    """
    Generate an operator file for the 4-state Hamiltonian of HBQ,
    suitable for the MCTDH code

    Input:
     modes           : a list of modes to be included in the Hamiltonian
     op_file         : the name of the operator file
     parameters_path : full path of the folder 'surface_parameters'
    """
    from os.path import join
    n_modes = len(modes)
    _eps    = 1e-6
    ind_0   = modes.index(0)
    
    if parameters_path is None:
        parameters_path = parameters_dir
    
    with open(op_file, 'w') as fh:
        
        # Header
        header = """
        OP_DEFINE-SECTION
        title
        4-state model (ground state and 3 coupled excited states) for HBQ
        end-title
        end-op_define-section
        
        """
        header = '\n'.join(line.lstrip() for line in header.splitlines())
        fh.write(header)
        fh.write('HAMILTONIAN-SECTION \n')
        fh.write('----------------------------------------------------\n')
        fh.write('modes | el ')
        for i in range(n_modes):
            modes_string = (i % 10 == 0)
            if modes_string: fh.write('\nmodes')
            fh.write(f' |  Q{(modes[i] + 1):<2} ')
        fh.write('\n----------------------------------------------------\n')
            
        # Kinetic energy operator
        for i in range(n_modes):
            m = modes[i]
            fh.write('%14.10f   |%-2i dq^2 \n' % (-0.5 * Omega[m,m], i + 2))
            for j in range(i + 1, n_modes):
                n = modes[j]
                if abs(Omega[m,n]) > _eps:
                   fh.write('%14.10f   |%-2i dq   |%-2i dq   \n' % (-Omega[m,n], i + 2, j + 2))
                   
        fh.write('\n')
        labels = ''
        
        # Diagonal Potentials        
        for s in range(4):
            # constant
            if 0 in modes:
                func_str = 'v' + str(s)
                fh.write('  1.0            |1 S%i&%i  |%-2i  %s \n' % (s + 1, s + 1, ind_0 + 2, func_str))
                label = func_str + ' = external1d{' + join(parameters_path, 'V' + str(s) + '.dat') + '}'
                labels += label + '\n'
            else:
                fh.write('%14.10f   |1 S%i&%i \n' % (V1D[s](0), s + 1, s + 1))
            
            # linear and quadratic terms
            for i in range(n_modes):
                n = modes[i]
                #
                if n == 0: continue
                #
                if 0 in modes:
                    # linear term
                    func_str = 'v' + str(s) + '_' + str(n + 1)
                    fh.write('  1.0            |1 S%i&%i  |%-2i  %s  |%-2i  q \n' % (s + 1, s + 1, ind_0 + 2, func_str, i + 2))
                    label = func_str + ' = external1d{' + join(parameters_path, 'Grad_' + str(s) + '_' + str(n + 1) + '.dat') + '}'
                    labels += label + '\n'
                    # quadratic term
                    func_str = 'v' + str(s) + '_' + str(n + 1) + '_' + str(n + 1)
                    fh.write('  0.5            |1 S%i&%i  |%-2i  %s  |%-2i q^2 \n' % (s + 1, s + 1, ind_0 + 2, func_str, i + 2))
                    label = func_str + ' = external1d{' + join(parameters_path, 'Hess_' + str(s) + '_' + str(n + 1) + '_' + str(n + 1) + '.dat') + '}'
                    labels += label + '\n'
                else:
                    fh.write('%14.10f   |1 S%i&%i |%-2i  q \n' % (grad[s][n + 1](0), s + 1, s + 1, i + 2))
                    fh.write('%14.10f   |1 S%i&%i |%-2i  q^2 \n' % (0.5 * Hess[s][n + 1,n + 1](0), s + 1, s + 1, i + 2))                
                # bilinear terms
                for j in range(i + 1, n_modes):
                    m = modes[j]
                    #
                    if m == 0: continue
                    #
                    if 0 in modes:
                        func_str = 'v' + str(s) + '_' + str(n + 1) + '_' + str(m + 1)
                        fh.write('  1.0            |1 S%i&%i  |%-2i  %s  |%-2i q  |%-2i q \n' % (s + 1, s + 1, ind_0 + 2, func_str, i + 2, j + 2))
                        label = func_str + ' = external1d{' + join(parameters_path, 'Hess_' + str(s) + '_' + str(n + 1) + '_' + str(m + 1) + '.dat') + '}'
                        labels += label + '\n'
                    else:
                        fh.write('%14.10f   |1 S%i&%i |%-2i  q |%-2i q \n' % (Hess[s][n + 1, m + 1], s + 1, s + 1, i + 2, j + 2))
                        
            # anharmonic terms (only if the mode n. 8 is included)
            if 8 in modes:
                ind_8 = modes.index(8)
                if 0 in modes:
                    func_str = 'v' + str(s) + '_9_9_9'
                    fh.write('  1.0            |1 S%i&%i  |%-2i  %s  |%-2i q^3 \n' % (s + 1, s + 1, ind_0 + 2, func_str, ind_8 + 2))
                    label = func_str + ' = external1d{' + join(parameters_path, 'c3_' + str(s) + '.dat') + '}'
                    labels += label + '\n'
                    #
                    func_str = 'v' + str(s) + '_9_9_9_9'
                    fh.write('  1.0            |1 S%i&%i  |%-2i  %s  |%-2i q^4 \n' % (s + 1, s + 1, ind_0 + 2, func_str, ind_8 + 2))
                    label = func_str + ' = external1d{' + join(parameters_path, 'c4_' + str(s) + '.dat') + '}'
                    labels += label + '\n'
                else:
                    fh.write('%14.10f   |1 S%i&%i |%-2i  q^3 \n' % (c3[s](0), s + 1, s + 1, ind_8 + 2))
                    fh.write('%14.10f   |1 S%i&%i |%-2i  q^4 \n' % (c4[s](0), s + 1, s + 1, ind_8 + 2))
                    
        # Inter-state couplings
        fh.write('\n')
        # 1-2 coupling
        if 0 in modes:
            func_str = 'v12'
            fh.write('  1.0            |1 S1&2  |%-2i  %s \n' % (ind_0 + 2, func_str))
            label = func_str + ' = external1d{' + join(parameters_path, 'coup12.dat') + '}'
            labels += label + '\n'
        else:
            fh.write('%14.10f   |1 S1&2 \n' % W12(0))
        # 1-3 and 2-3 couplings
        for i_coup in [1, 2, 3, 4]:
            if i_coup in modes:
                ind_coup = modes.index(i_coup)
                if 0 in modes:
                    func_str = 'v13_' + str(i_coup + 1)
                    fh.write('  1.0            |1 S1&3  |%-2i  %s  |%-2i  q \n' % (ind_0 + 2, func_str, ind_coup + 2))
                    label = func_str + ' = external1d{' + join(parameters_path, 'coup13_' + str(i_coup) + '.dat') + '}'
                    labels += label + '\n'
                    #
                    func_str = 'v23_' + str(i_coup + 1)
                    fh.write('  1.0            |1 S2&3  |%-2i  %s  |%-2i  q \n' % (ind_0 + 2, func_str, ind_coup + 2))
                    label = func_str + ' = external1d{' + join(parameters_path, 'coup23_' + str(i_coup) + '.dat') + '}'
                    labels += label + '\n'
                else:
                    fh.write('%14.10f   |1 S1&3 |%-2i  q \n' % (W13[i_coup - 1](0), ind_coup + 2))
                    fh.write('%14.10f   |1 S2&3 |%-2i  q \n' % (W23[i_coup - 1](0), ind_coup + 2))
            
        ##########
        fh.write('end-hamiltonian-section \n\n')
        fh.write('LABELS-section \n')        
        fh.write(labels)
        fh.write('end-labels-section')
        fh.write('\n\n')
        fh.write('end-operator')
               

        
        
        
    
    
