
import numpy as np
import qutip
import os
import pandas as pd

#Negativity Helpers

def negativity(density_data):
    """
    Calculates the negativity over a splitting partition for the given density matrix data
    NOTE: assumes ring topology of Ising Model
    
    Parameters
    ----------
        density_data : array_like
            2^(Np) x 2^(Np) numpy array that stores the data for the density matrix for the primary system

    """
    N = int(np.log2(len(density_data)))

    rho = qutip.Qobj(inpt = density_data, dims=[[2]*N, [2]*N])
    splitting_subspace = ([1]* (N//2)) +  ([0] * (N//2))

    if N %2 == 1:
        splitting_subspace += 1

    eigenvalues = qutip.partial_transpose(rho, mask = splitting_subspace).eigenenergies()
    
    return __calc_negativity(eigenvalues)

 

def __calc_negativity(eig):
    """
    Sums the negative eigenvalues
    """
    neg = 0
    for val in np.real(eig):
        if val < 0:
            neg += abs(val)

    return neg