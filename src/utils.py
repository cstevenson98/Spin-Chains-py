import itertools
import sys
import scipy.sparse as SS
import numpy as np



# # FUNCTION: Nest_flatten
    
# Take a structured list of sub-matrices and remove their 'site label'.

def Nest_flatten(nested):
    out = []
    for S in nested:
        site = []
        for block in S:
            for eig_op in block:
                site.append(tuple(eig_op))
        out.append(site)
                
    return out


# # FUNCTION: dm_shape

# DESCRIPTION: Tensor structure for a lattice of N 2LS.

def dm_shape(N):
    T = [2 for i in range(N)]
    return [T, T]
    
# # FUNCTION: super_shape

# DESCRIPTION: The correct tensor structure for the open system consisting of N 2LS.

def super_shape(N):
    T = [2 for i in range(N)]
    return [[T, T], [T, T]]


# FUNCTION: block_to_2N

# DESCRIPTION: Takes a matrix lying on the mth diagonal of a block
#              matrix and converts it to a 2^N x 2^N sparse matrix.
    
def block_to_2N(block, N, x_pos, y_pos):

    out = SS.dok_matrix((2**N, 2**N))
    for i in range(len(block)):
        for j in range(len(block[0])):
            out[x_pos+i, y_pos+j] = block[i][j]
            
    return out

# FUNCTION: block_visualisation

def block_visualisation(block):
    for i in range(len(block)):
        for j in range(len(block[0])):
            if block[i][j] == 0.:
                sys.stdout.write(' ')
            else:
                sys.stdout.write('#')
        sys.stdout.write('\n')
    sys.stdout.write('\n')


def restart_line():
    sys.stdout.write('\r')
    sys.stdout.flush()

def loading_bar(k, toolbar_width):
    sys.stdout.write("["+int(20*(k/toolbar_width))*"="+(20-int(20*(k/toolbar_width)))*" "+"]")
    sys.stdout.write("  {0:.0f} %".format(float(k)/toolbar_width * 100))
    sys.stdout.write("  ({0}/{1})".format(k, toolbar_width))
    sys.stdout.flush()
    restart_line()
