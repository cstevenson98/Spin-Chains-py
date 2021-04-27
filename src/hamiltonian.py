import numpy as np

from utils import *
from spin_alg import *
from block_matrices import *

# # FUNCTION: HXY

# DESCRIPTION: Derive the pure XY spin-chain Hamiltonian, H, with no
#              external field.
#              H = Sum[s^+_i*s^-_(i+1)+h.c.] + g*Sum[s^z_i].

def HXY(N, g, boundary_conditions):
    out = []
    for r in range(N+1):
        basis_r = Kbits(N, r)
        dim_r = len(basis_r)
        H_r = np.zeros((dim_r, dim_r),dtype=float)
        
        mapped = list(map(lambda x: Hop(x, boundary_conditions), basis_r))

        for k in range(dim_r):
            if mapped[k]:
                for result in mapped[k]:
                    H_r[k, basis_r.index(result)] += 1
        
        for i in range(dim_r):
            H_r[i,i] += g*(2*r-N)
            
        out.append(H_r)

    return out

def HXY_block(N, g, boundary_conditions):
    out = Block_Matrix(N, 0)
    out.data = HXY(N, g, boundary_conditions)
    return out
