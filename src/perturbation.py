import itertools

import scipy.sparse as SS
import scipy.misc as SM
import numpy as np
import numpy.linalg as npLA

from qutip import *
from utils import *
from hamiltonian import *
from block_matrices import *

# # FUNCTION: perturbationV

# DESCRIPTION: The operator which is the pair creation/annihilation
#              operator sum_{i=1}^N {sig^+_isig^+_{i+1}} or H.c.,
#              depending on input.

def perturbationV(kind, N):
    S_minus_minus_list = []
    for i in range(N-1):
        S_minus_minus_list.append(S_minus_minus_block(N, i, i+1))
    S_minus_minus_list.append(S_minus_minus_block(N, N-1, 0))

    if kind == "--":
        return sum(S_minus_minus_list)
    elif kind == "++":
        return sum(map(lambda x: x.transpose(), S_minus_minus_list))
    

# # FUNCTION: corrections

# DESCRIPTION: Given a set of pertubative operators and the
#              eigenspectrum (both eigenenergies and eigenstates),
#              calculate the set of perturbative corrections which go
#              into the wavefunction.

def corrections(Vset, e_en, e_vec):
    # For each of the operators in Vset, calculate the block matrix of
    # matrix elements of the eigenvectors.
    correct_set = []
    for V in Vset:
        correct_set.append(e_vec*V*e_vec.transpose())

    # Now calculate the full corrections by dividing by the energy
    # difference.
    for corr in correct_set:
        corr.shape[1]
        for (i, block) in enumerate(corr.data):
            for ((m, n), val) in np.ndenumerate(block):
                if (val != 0):
                    if corr.shape[1] >= 0:
                        corr.data[i][m][n] = val/(e_en[i][m] -
                                                  e_en[i+corr.shape[1]][n])
                    elif corr.shape[1] < 0:
                        corr.data[i][m][n] = val/(e_en[i-corr.shape[1]][m]
                                                  - e_en[i][n])

    return correct_set
