import sys

import timing

from qutip import *
from spin_alg import *
from scipy.sparse import hstack
import numpy as np
import numpy.linalg as LA
import itertools

N = int(sys.argv[1])
J=1.
deltax=1.
g=float(sys.argv[2])

# Some functions to generate configuration basis.
def binary_str(n):
    return  list(itertools.product([0,1], repeat = n))

def num_of_ones(bin_str):
    out = []
    for i in bin_str:
        num_ones=0
        for j in i:
            num_ones += j
        out.append(num_ones)
    return out

def config_basis(tup):
    ket = []
    for i in tup:
        if i==1:
            ket.append(basis(2,1))
        else:
            ket.append(basis(2))
    return tensor(ket)

bin_states = binary_str(N)
r_ordered_states = []
for i in np.argsort(num_of_ones(bin_states)):
    r_ordered_states.append(config_basis(bin_states[i]))


sx_list = []
sz_list = []
sm_list = []
sp_list = []
for i in range(N):
    oper_list = []
    for j in range(N):
        oper_list.append(identity(2))
    oper_list[i] = sigmax()
    sx_list.append(tensor(oper_list))
    oper_list[i] = sigmaz()
    sz_list.append(tensor(oper_list))
    oper_list[i] = sigmam()
    sm_list.append(tensor(oper_list))
    oper_list[i] = sigmap()
    sp_list.append(tensor(oper_list))

H=0.
for i in range(N-1):
    H += J*(sm_list[i]*sp_list[i+1] + sp_list[i]*sm_list[i+1]) + g*sz_list[i]
H += g*sz_list[N-1]
H += J*sp_list[N-1]*sm_list[0] + J*sm_list[N-1]*sp_list[0]

# Definition of Liouvillian
L_check = sum([2*lindblad_dissipator(sm_list[i]) for i in range(N)])

T = basis_mat(N)
H_check = LA.inv(T).dot((H.data.todense())).dot(T)

def init(m,i):
    if m == i:
        return basis(2,0)
    else:
        return basis(2,1)

def expec(m,i):
    if m == i:
        return sigmaz()
    else:
        return identity(2)
    
def dm_shape(m):
    A = [2 for i in range(m)]
    return [A, A]
    
# times = np.linspace(0.0, 10.0, 100)
# psi0 = rand_dm(2**N,dims=dm_shape(N))

# result = mesolve(L_check, psi0, times, [], [tensor([expec(0, i) for i in range(N)])]) 
# fig, ax = subplots()
# ax.plot(result.times, result.expect[0])
# ax.plot(result.times, result.expect[0])
# ax.set_xlim([0,2])
# ax.set_ylim([-1.2,1.2])
# show()

print(L_check)
