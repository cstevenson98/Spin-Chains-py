#!/home/conor/anaconda3/bin/python

import sys
import getopt
import itertools
import scipy.sparse as SS
import numpy as np
import numpy.linalg as npLA
import multiprocessing as mp

import timing

from qutip import *
from utils import *
from hamiltonian import *
from master import *
from block_matrices import *
from perturbation import *

###############################################################
# PARSING COMMAND LINE ARGUMENTS (PERTURBATION THEORY OR NOT) #
###############################################################
#
# WAYS TO USE THIS PROGRAM:
#
# --- FULLY UNPERTURBED: ---
#
#     python main.py -a N,g,omega0,b 
#
# --- WITH PERTURBATION THEORY: ---
#
#     python main.py -a N,g,omega0,b -p Lambda
#
#
# --- UNPERTURBED DISSIPATOR, PERTURBED UNITARY EVOLUTION: ---
#
#     python main.py -a N,g,omega0,b --upd Lambda
#
#     NB: This is the Liouvillian obtained by calculating the
#     dissipation using the eigenstates of the unperturbed Hamiltonian
#     H_0, but while using evolving unitarily with H_0 + H_I.
#
# --- OPTIONS:
#
#     --multiprocessing=no_of_cores
#
#        This will split the calculation of the Liouvillian across
#        multiple no of cores.
#
###############################################################

try:
    opts, args = getopt.getopt(sys.argv[1:], "ha:p:u:", ["upd=", "multiprocessing=", "PBC"])
except getopt.GetoptError:
    print('main.py -a <N,g,omega0,b> [-pu <lambda>]')
    sys.exit(2)

flag = 'Unperturbed'
no_of_cores = 1
boundary_conditions = 'OBC'

for opt, arg in opts:
    if opt == '-h':
        print('main.py -a <N,g,omega0,b> [-pu <lambda>]')
        sys.exit()
    elif opt == '-a':
        arg = arg.split(",")
        N = int(arg[0])
        g = float(arg[1])
        omega0 = float(arg[2])
        b = float(arg[3])
    elif opt == '-p':
        flag = 'PerturbationTheory'
        Lambda = float(arg)
    elif opt in ("-u", "--upd"):
        flag = 'UnperturbedDissipator'
        Lambda = float(arg)
    elif opt == "--multiprocessing":
        no_of_cores = int(arg)
    elif opt == "--PBC":
        boundary_conditions = 'PBC'

        
####################
# MAIN CALCULATION #
####################

def Spec_func(w, w0, b):
    return np.exp(-((w-w0)/b)**2)

J = lambda x: Spec_func(x, omega0*bnd_wid, b*bnd_wid)

T = basis_mat(N)
T_i = npLA.inv(T)

def liou(x):
    return master(x, J, N)

if __name__ == "__main__":
    
    # Begin by constructing the hamiltonian using binary string manip.
    sys.stdout.write('Constructing Hamiltonian...')
    H = HXY_block(N, g, boundary_conditions)
    sys.stdout.write(' DONE.\n')
    print(H)

    
    # Since number conserving, block diagonal structure. So can simply
    # map over nested struct. Specifically need to use np.linalg.eigh
    # for hermitian operator (guarantees orthonormality).
    sys.stdout.write('Calculating eigenstates...')
    eig_block_list = list(map(npLA.eigh, H.data))
    sys.stdout.write(' DONE.\n')

    energies = [eig_block_list[i][0] for i in range(N+1)]
    en = np.concatenate(np.array(energies))

    # Eigenvectors are stored in the rows, we turn them into columns
    # in order to represent kets.
    e_vectors = [eig_block_list[i][1].T for i in range(N+1)]
    e_vectors = Block_Matrix(N, 0, mydata = e_vectors)

    # # Calculate bandwidth: The sweep of a gaussian spectral function
    # # should make sure to cover the range of transitions allowed by the
    # # terms in a master equation, so omega0 is in units of the bandwidth
    # # which we now calculate:
    freq_diff = map(lambda x: x[1]-x[0], itertools.product(en, en))
    bnd_wid = max(freq_diff)

    # The spin lower operators are formed similarly to the
    # Hamiltonian, by manipulating binary strings and using .key to
    # encode the matrix form in the ordered spin basis (we will need
    # to transform into the basis QuTIP uses later.
    sys.stdout.write('Forming spin-lowering operators...')
    S_list = [S_minus_block(N, i) for i in range(N)]
    sys.stdout.write(' DONE.\n')


    # The filename suffix which containes the system parameters.
    filename = "N="+str(N)+"-g="+str(g)+"-omega0="+str(omega0)+"-b="+str(b)

    if flag == 'PerturbationTheory':

        # In the perturbation approach, we construct first the
        # corrections, then the eigenoperators formed in the perturbed
        # basis. Then it is a matter of pairing these eigenoperators
        # up correctly as specified by the Born Markov master
        # equation. This job can be split amongst N cores, as once all
        # of the Eigenoperators are calculated, on simply needs to
        # promote them to super-operators, so since this is the most
        # computationally time intensive (I think) then it is useful
        # at this stage to parallelarise the loop.

        # Fist calculate the perturbations due to a pair of spin
        # raising and spin lowering operators.
        S_plusplus = perturbationV("++", N)
        S_minuminu = perturbationV("--", N)

        # This allows one to then calculate the perturbative
        # corrections in matrix form. (It is at this stage that one
        # needs to be careful about degeneracies. These can be avoided
        # if |g| is sufficiently large, since in that limit the bands
        # are well seperated, and will not overlap at all. Since the
        # perturbations are pairs of raising and lowering, we never
        # get energy differences from within the same
        # band. Additionally, the perturbation series is most accurate
        # here as it is well within the perturbative assumption (g >>
        # J, or something.)
        sys.stdout.write('Calculating perturbative corrections...')
        corr = corrections([S_plusplus, S_minuminu], energies, e_vectors)
        sys.stdout.write(' DONE.\n')

        # The eigenoperators are all calculated here, using the
        # complicated formulae derived in the note. There are a set of
        # eigenoperators for each site in the spin-chain.
        sys.stdout.write('Calculating perturbative eigenoperators...')
        e_ops = pert_eigenoperators(S_list, energies, e_vectors, Lambda, corr)
        sys.stdout.write(' DONE.\n')

    elif flag == 'Unperturbed' or flag == 'UnperturbedDissipator':
        # Do non-perturbed calculation
        sys.stdout.write('Calculating eigenoperators of H...')
        e_ops = eigenoperators(S_list, energies, e_vectors)
        sys.stdout.write(' DONE.\n')
        
    
    # Now we convert all the eigenoperators to QuTiP quantum
    # objects, in the QuTiP basis. We create a copy of the
    # eigenoperators multiplied by the spectral function evaluated
    # at the energy difference of the eigenoperators.
    Jeops = [map(lambda x: Qobj(T.dot((J(x[0])*x[1]).dot(T_i))), e_ops[p])
             for p in range(len(e_ops))]
    Eops = [map(lambda x: Qobj(T.dot((x[1]).dot(T_i))), e_ops[p])
            for p in range(len(e_ops))]

    # Within the set of eigenoperators on each site, create all pairs.
    e_ops = [list(itertools.product(Eops[p], Jeops[p])) for p in range(len(e_ops))]

    # We can now flatten this list as all that remains is to
    # iterate through all of the pairs and 
    flat_e_ops = [pair for site in e_ops for pair in site]

    # Split up each set of pairs amonst 'no_of_cores' processors
    processes = mp.Pool(no_of_cores)
    e_ops_multiprocessing = [flat_e_ops[p::no_of_cores] for p in range(no_of_cores)]

    # Finally we call master for each of these 
    if flag == "PerturbationTheory":
        sys.stdout.write('Calculating perturbed Liouvillian...')
    elif flag == "Unperturbed" or flag == "UnperturbedDissipator":
        sys.stdout.write('Calculating Liouvillian...')

    L = sum(processes.map(liou, e_ops_multiprocessing))
    sys.stdout.write(' DONE.\n')
    L.dims = super_shape(N)

    if flag == "PerturbationTheory":
        filenameL = "L::PERT::" + filename
        filenameD = "D::PERT::" + filename
        filenamerhoSS = "rhoSS::PERT::" + filename
    elif flag == 'Unperturbed':
        filenameL = "L" + filename
        filenameD = "D" + filename
        filenamerhoSS = "rhoSS" + filename
    elif flag == 'UnperturbedDissipator':
        filenameL = "L::U-PDiss::" + filename
        filenameD = "D::U-PDiss::" + filename
        filenamerhoSS = "rhoSS::U-PDiss::" + filename

    

    # Let us finally put the Hamiltonian of the system, which is
    # either perturbed or unperturbed in the 
    if flag == 'PerturbationTheory' or flag == 'UnperturbedDissipator':
        fullH = H.to_sparse() + Lambda*(perturbationV("++", N).to_sparse() +
                                        perturbationV('--', N).to_sparse())
    else:
        fullH = H.to_sparse()
    fullH = Qobj(T.dot(fullH.dot(T_i)))
    
    fullH.dims = dm_shape(N)

    #Finally we add the unitary evolution to the master equation and
    #solve for the steady state.
    
    sys.stdout.write("Calculating steady state solution... ")
    
    i_H_comm_rho = -liouvillian(fullH)
    rho_dot = L-i_H_comm_rho
    rhoSS = steadystate(rho_dot)
    sys.stdout.write("DONE.\n")

    # Save these quantum objects -- much analysis can be done with these files.
    qsave(rho_dot, filenameL)
    qsave(L, filenameD)
    qsave(rhoSS, filenamerhoSS)

    # Exit newline
    sys.stdout.write('\n')
