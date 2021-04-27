# # # MODULE: spin_alg
# # # # # # # # # # # # # # # # #
# # DESCRIPTION: This module aims to calculate spin operators in matrix form,
# #              and calculate, for example, the Hamiltonian of the XY model.
# #
# # Conor Stevenson 2017
#

import itertools
import numpy as np
import scipy as sp
import scipy.linalg as LA
import scipy.misc as SM
import scipy.sparse as SS

from block_matrices import *

#from qutip import shape

# # FUNCTION: Hop

# DESCRIPTION: Determines all strings resulting from exactly one
# "01"<->"10" flip. E.g.: hop('101')=['011', '110'],
# hop('100') = ['010', '001'].

def Hop(a, boundary_conditions):
    N = len(a)
    out = []
    for i in range(N):
        temp = list(a)
        if not(i == N-1):
            if not(a[i] == a[i+1]):
                store = a[i]
                temp[i] = temp[i+1]
                temp[i+1] = store
                out.append(''.join(temp))
        if i == N-1 and boundary_conditions == 'PBC':
            if not(a[N-1] == a[0]):
                store = a[N-1]
                temp[N-1] = temp[0]
                temp[0] = store
                out.append(''.join(temp))

    return out


# # FUNCTION: Kbits

# DESCRIPTION: Generate all binary strings of length n with exactly k 1s.

def Kbits(N, r):
    out = []
    for bits in itertools.combinations(list(range(N)), r):
        s = ['0'] * N
        for b in bits:
            s[b] = '1'
        out.append(''.join(s))
    return out


# # FUNCTION: Down

# DESCRIPTION: If position l is '1', flips to '0', else returns nothing.

def Down(bstr, l):
    out = list(bstr)
    if not out[l] == '1':
        return 0
    else:
        out[l] = '0'
    return ''.join(out)


# # FUNCTION: S_minus

# DESCRIPTION: Contruct the N off-diagonal block matrices which define
# the spin-lowering operator in the ordered basis.

def S_minus(size, l):
    N = size
    
    if l < 0 or l >= N:
        print('S_MINUS: Error, incorrect site label')
        return 0
              
    out = []
    for k in range(N):
        basis_k, basis_k1 = Kbits(N, k), Kbits(N, k+1)

        dim1 = len(basis_k)
        dim2 = len(basis_k1)
        sig = np.zeros((dim1, dim2),dtype=int)
        
        count = 0
        for col in basis_k1:
            flipped = Down(col, l)
            if not flipped == 0:
                sig[basis_k.index(flipped), count] += 1
            count += 1
        out.append(sig)
        
    return out


# # FUNCTION: S_minus_minus

# DESCRIPTION: The (N-1) blocks of the sig^-_l sig^-_m

def S_minus_minus(size, l, m):
    N = size

    if l < 0 or l >= N:
        print('S_MINUS_MINUS: Error, incorrect site label')
        return 0

    out = []
    for k in range(N-1):
        basis_k, basis_k2 = Kbits(N, k), Kbits(N, k+2)

        dim1 = len(basis_k)
        dim2 = len(basis_k2)
        sig = np.zeros((dim1, dim2), dtype=int)

        count = 0
        for col in basis_k2:
            flipped = Down(col, m)
            if not flipped == 0:
                flipped = Down(flipped, l)
                if not flipped == 0:
                    sig[basis_k.index(flipped), count] += 1
            count += 1
        
        out.append(sig)

    return out


# # FUNCTION: S_minus_block

# DESCRIPTION: Creates a block matrix object describing spin lowering on one site

def S_minus_block(size, l):
    out = Block_Matrix(size, 1)
    out.data = S_minus(size, l)
    return out

def S_minus_minus_block(size, l, m):
    out = Block_Matrix(size, 2)
    out.data = S_minus_minus(size, l, m)
    return out


## FUNCTION: basis_mat

# DESCRIPTION: The basis transformation matrix from the spin-ordered basis
#              (e.g.:{'000','100','010','001','110','101','011','111'}) to
#              the tensor basis {ijk : Binary number ijk decreases}.

def basis_mat(N):
    sorted_basis = []
    for i in range(N+1):
        for s in Kbits(N,i):
            sorted_basis.append(s)

    bin_down = []
    for byte in itertools.product(['1','0'], repeat = N):
        out = ''
        for bit in byte:
            out = out + bit
        bin_down.append(out)
            
    keys = []
    for string in bin_down:
        keys.append(sorted_basis.index(string))

    def delt(i,j):
        if i==j:
            return 1
        else:
            return 0

    idd = [[delt(i,j) for i in range(2**N)] for j in range(2**N)]


    basis_change = []
    for k in keys:
        basis_change.append(idd[k])


    return np.array(basis_change)

    
