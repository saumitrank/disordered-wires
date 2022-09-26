"""
This module contains utility functions
"""

import matplotlib.pyplot as plt
import numpy as np
import kwant




def plot_band_structure(lead:kwant.system.InfiniteSystem, momenta:np.array)->tuple:
    """Plots band structure given lead and momenta

    Args:
        lead (kwant.system.InfiniteSystem): Lead defined in KWANT
        momenta (np.array): Points to find band structure

    Returns:
        tuple: figure and matplotlib axis
    """

    #Get intra and inter unit-cell Hamiltoninan
    ham = lead.cell_hamiltonian()
    hop = np.empty(ham.shape, dtype=complex)
    t = lead.inter_cell_hopping()
    hop[:, : t.shape[1]] = t
    hop[:, t.shape[1] : ] = 0

    momenta = np.array(momenta)
        
    #Define full Hamiltonian in k-space
    H = ham[np.newaxis, :] \
        + hop[np.newaxis, :]*np.exp(-1j*momenta)[:, np.newaxis, np.newaxis] \
        + hop.T.conjugate()[np.newaxis, :]*np.exp(1j*momenta)[:, np.newaxis, np.newaxis]
    
    energies = np.linalg.eigvalsh(H)
    
    #Plot
    fig, ax = plt.subplots()
    ax.plot(momenta, energies)

    ax.set_xlabel(r"$k_x \;  \; [a^{-1}]$")
    ax.set_ylabel(r"$E(k_x) \; \; [t]$")
    ax.set_xlim(np.min(momenta), np.max(momenta))
    
    return fig, ax