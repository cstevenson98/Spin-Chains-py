import numpy as np
import scipy.sparse as SS

def fact(n):
    if n==0: return 1
    return n*fact(n-1)

def bin(n, r):
    return fact(n)//(fact(n-r)*fact(r))


# # CLASS: Block_Matrix

# DESCRIPTION: A class which defines a binomial block matrix. This is
#              a 2^N x 2^N matrix, with blocks which range through the
#              binomial coefficients, (N C r). We may have these
#              blocks offset from the main diagonal by an amount
#              denoted m. So if the blocks are on the main diagonal,
#              m=0. Then there are N+1 blocks, with shapes (N C 0)x(N
#              C 0), (N C 1)x(N C 1),... (N C (N-1))x(N C (N-1)), (N C
#              N)x(N C N).
#
#              If m is non-zero, then the blocks are not square and
#              lie off the main diagonal, with shapes (N C 0)x(N C m),
#              (N C 1)x(N C m+1), ... (N C N-m)x(N C N).

class Block_Matrix:

    # METHOD: __init__
    
    # DESCRIPTION: Create an empty binomial block matrix, offset by
    #              'offset' from the main diagonal, with dimension
    #              2^size x 2^size.
    
    def __init__(self, size, offset, mytype = int, mydata = 0):

        if offset > size or offset < -size:
            raise Exception('Invalid main-diagonal offset.')

        self.shape = [size, offset]

        if mydata == 0:
            if offset >= 0:
                self.data = [np.zeros((bin(size, i), bin(size, offset+i)),
                                      dtype=mytype) for i in range(size-offset+1)]
            else:
                self.data = [np.zeros((bin(size, -offset+i), bin(size, i)),
                                      dtype=mytype) for i in range(size+offset+1)]
        else:

            # Check if handed data is of the correct type. I.e. check
            # that is a list of lists, and that the sub-lists have the
            # correct shapes for a binomial-block-matrix.
            if offset >= 0:
                shapes = [(bin(size, i), bin(size, offset+i))
                          for i in range(size-offset+1)]
            else:
                shapes = [(bin(size, -offset+i), bin(size, i))
                          for i in range(size+offset+1)]
                
            for i in range(size-abs(offset)+1):
                if not(type(mydata[i]) is np.ndarray):
                    raise Exception(
                        "Value error: The input data is not a list of Numpy arrays")
                if not(np.shape(mydata[i]) == shapes[i]):
                    raise Exception(
                        "Value error: The input data is not of the correct form")
            self.data = mydata


    # METHOD: __str__
            
    # DESCRIPTION: Method to tell python how to print the object.

    def __str__(self):
        out = ''
        out += 'Binomial block matrix object:\n'
        out += 'Shape = ' + str(self.shape) + '\n'
        out += 'Data = \n'
        for i in range(self.shape[0]-abs(self.shape[1])+1):
            out+=str(self.data[i])+'\n'
        return out


    # METHOD: __repr__
            
    # DESCRIPTION: Method to tell python how to print the object in an
    #              interactive setting.

    def __repr__(self):
        out = ''
        out += 'Binomial block matrix object:\n'
        out += 'Shape = ' + str(self.shape) + '\n'
        out += 'Data = \n'
        for i in range(self.shape[0]-abs(self.shape[1])+1):
            out+=str(self.data[i])+'\n'
        return out


    # METHOD: __add__

    # DESCRIPTION: Overload the addition operator. By convension, we
    #              do not allow addition of differently shaped block
    #              structures, as this would result in a hybrid
    #              structure.

    def __add__(self, other):
        if other == 0:
            return self
        
        if self.shape != other.shape:
            raise Exception('Block matrix shape mismatch.')
        
        out = self
        for i in range(self.shape[0]-abs(self.shape[1])+1):
            out.data[i] += other.data[i]
            
        return out

    
    # METHOD: __radd__

    # DESCRIPTION: Reverse addition. Allows implementation of, for
    #              example the 'sum' function for block matrices.
    
    def __radd__(self, other):
        return self + other

    
    # METHOD: __sub__

    # DESCRIPTION: Overload the addition operator. By convension, we
    #              do not allow addition of differently shaped block
    #              structures, as this would result in a hybrid
    #              structure.

    def __sub__(self, other):
        if self.shape != other.shape:
            raise Exception('Block matrix shape mismatch.')
        
        out = self
        for i in range(self.shape[0]-abs(self.shape[1])+1):
            out.data[i] -= other.data[i]
            
        return out


    # METHOD: __rmul__

    # DESCRIPTION: Define multiplication of matrices by numbers.

    def __rmul__(self, number):
        out = Block_Matrix(self.shape[0], self.shape[1], mytype = type(number))

        for i in range(self.shape[0]-abs(self.shape[1])+1):
            out.data[i] = number*self.data[i]

        return out

    
    # METHOD: __mul__

    # DECRIPTION: Overload of the multiplication function for block
    #             matrices.
    #           
    #             If the object on the left is not a block structure,
    #             then revert back to the __rmul__ routine, which
    #             defines multiplication by numbers.

    def __mul__(self, other):
        if not(type(other) is Block_Matrix):
            return other*self

        if self.shape[0] != other.shape[0]:
            raise Exception('Block matrix shape mismatch.')

        size = self.shape[0]
        M = self.shape[1]
        N = other.shape[1]
        MplusN = self.shape[1]+other.shape[1]
        
        out = Block_Matrix(self.shape[0], MplusN)

        # Split into cases: m >= n (m & n >= 0; m > 0 & n < 0; etc.) and
        #                   m < n (sim. )
        if abs(MplusN) > out.shape[0]:
            out.shape = (0, 0)
            out.data = 0
            return out

        if M >= 0:
            if N >= 0:
                for i in range(len(out.data)):
                    out.data[i] = np.dot(self.data[i], other.data[M+i])

        if M >= 0:
            if N < 0:
                if abs(M) >= abs(N):
                    for i in range(len(out.data)-abs(N)):
                        out.data[i] = np.dot(self.data[i],
                                             other.data[abs(M)-abs(N)+i])
                
                if abs(M) < abs(N):
                    for i in range(len(out.data)-abs(M)):
                        out.data[i] = np.dot(self.data[abs(N)-abs(M)+i],
                                             other.data[i])

        if M < 0:
            if N >= 0:
                if abs(N) >= abs(M):
                    for i in range(abs(M), len(out.data)):
                        out.data[i] = np.dot(self.data[i-abs(M)],
                                             other.data[i-abs(M)])

                if abs(N) < abs(M):
                    for i in range(abs(N), len(out.data)):
                        out.data[i] = np.dot(self.data[i-abs(N)],
                                             other.data[i-abs(N)])

        if M < 0:
            if N < 0:
                for i in range(len(out.data)):
                    out.data[i] = np.dot(self.data[-N+i], other.data[i])

        return out

    
    # METHOD: transpose

    # DESCRIPTION: Take the transpose of the block matrix.
        
    def transpose(self):
        out = Block_Matrix(self.shape[0], -self.shape[1])
        for i in range(self.shape[0]-abs(self.shape[1])+1):
            out.data[i] = np.transpose(self.data[i])
        return out


    # METHOD: to_sparse

    # DESCRIPTION: Returns a 2^size-dimensional square matrix
    #              containing the content of the block matrix.

    def to_sparse(self):

        N = self.shape[0]
        M = self.shape[1]

        N_C_k = [bin(N, k) for k in range(N+1)]
        N_C_k_Cummulative = [0]
        for k in range(1,N+1):
            N_C_k_Cummulative.append(sum(N_C_k[:k]))

        coords = []
        if self.shape[1] >= 0:
            for i, j in enumerate(zip(N_C_k_Cummulative[:N+1],
                                      N_C_k_Cummulative[self.shape[1]:])):
                    coords.append(j)
        else:
            for i, j in enumerate(zip(N_C_k_Cummulative[-self.shape[1]:],
                                      N_C_k_Cummulative[:N+1])):
                coords.append(j)

        out = SS.dok_matrix((2**N, 2**N))

        for k in range(len(self.data)):
            mat = self.data[k]
            for i in range(len(mat)):
                for j in range(len(mat[0])):
                    out[coords[k][0]+i, coords[k][1]+j] = mat[i][j]

        return out
