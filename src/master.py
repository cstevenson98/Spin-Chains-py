import itertools

import scipy.sparse as SS
import scipy.misc as SM
import numpy as np
import numpy.linalg as npLA
from itertools import *


from qutip import *
from utils import *
from hamiltonian import *
from block_matrices import *



################################################################
################## UNPERTURBED NUMERICS ########################
################################################################

# # FUNCTION: eigenoperators

# DESCRIPTION: Utilise the block matrix structures

def eigenoperators(opers, eig_en, eig_v):
    # H - Hamiltonian in block-matrix.
    # opers - List of block-matrices of operators to
    #         calculate eigenoperators for.
    # eig_en - Eigen-energies
    # eig_v - Eigen-vectors in block-matrix form

    # Calculate the matrix elements
    mat_el_list = []
    for op in opers:
        mat_el_list.append(eig_v*op*eig_v.transpose())

    eop = []
    for site in mat_el_list:
        eop_site = []
        for (i, block) in enumerate(site.data):
            eop_block = []
            for ((m, n), val) in np.ndenumerate(block):
                op = val*np.outer(eig_v.data[i][m], eig_v.data[i+1][n])
                omega = eig_en[i][m]-eig_en[i+1][n]
                eop_block.append((op, omega))

            eop_site.append(eop_block)

        eop.append(eop_site)

    N = len(eop)
    
    N_C_k = [bin(N, k) for k in range(N+1)]
    N_C_k_Cummulative = [0]
    for k in range(1,N+1):
        N_C_k_Cummulative.append(sum(N_C_k[:k]))

    # Join all eigenoperators on each site, in preparation for
    # sum in master equation. First loop over the sites, then the
    # blocks in each site. Convert the operators in each block into
    # a full 2^N square matrix and convert to QuTiP Qobj object.
    # Perform a basis transformation on these blocks also.
    new_eops = []
    for eops_site in eop:
        new_eops_site = []
        for (i, mat_omega_list) in enumerate(eops_site):
            x_pos = N_C_k_Cummulative[i]
            y_pos = N_C_k_Cummulative[i+1]
            for (j, (mat, omega)) in enumerate(mat_omega_list):
                block = block_to_2N(mat, N, x_pos, y_pos).todense()
                new_eops_site.append((omega, block))
        new_eops.append(new_eops_site)
        
    return new_eops

################################################################
################## PERTURBATION THEORY #########################
################################################################


# # FUNCTION: pert_eigenoperator

# DESCRIPTION: Calulates the perturbed eigenoperators and outer
#              products which go into the perturbative Born-Markov
#              Master equation.

def pert_eigenoperators(opers, eig_en, eig_v, Lambda, corrections):
    N = len(eig_en)-1

    # Extract the correction for plusplus and minusminus:
    PP = corrections[0]
    MM = corrections[1]
    pp = PP.data
    mm = MM.data

    # Calculate all of the matrix elements, which begin with index 0,
    # linking the 0-spin subspace to the 1-spin subspace, up to the
    # (N-1)th index, containing the matrix linking the (N-1)-spin
    # subspace to the N-spin subspace.

    mat_el_list = []
    for op in opers:
        mat_el_list.append(eig_v*op*eig_v.transpose())

    pert_mat_elems_list = []
    for mat_el_block in mat_el_list:
        L = len(mat_el_block.data)
        
        tuplemult = lambda x: x[0].dot(x[1])
        def zip_and_mult(A, B):
            zipped = zip(A, B)
            mult = map(tuplemult, zipped)
            return list(mult)

        # Perturbative corrections for offset of one (n = m-1)
        sigP = zip_and_mult(mat_el_block.data[1:], pp)
        Msig = zip_and_mult(map(lambda x: x.transpose(), mm), \
                            mat_el_block.data[:L-1])

        # Perturbative correction for offset of three (n = m+3)
        sigM = zip_and_mult(mat_el_block.data[:L-2], mm[1:])
        Psig = zip_and_mult(map(lambda x: x.transpose(), pp[:len(pp)-1]), \
                            mat_el_block.data[2:])
        
        tupleadd = lambda x: x[0] + x[1]
        def zip_and_add(A, B):
            zipped = zip(A, B)
            added = map(tupleadd, zipped)
            return list(added)

        # Perturbative matrix elements of sigma minus for n = m-1        
        sigP_plus_Msig_zipped = zip_and_add(sigP[1:], Msig[:len(Msig)-1])
        Pert_mat_el_offsetbyone = [sigP[0]] + sigP_plus_Msig_zipped + \
                                  [Msig[len(Msig)-1]]

        # Perturbative matrix elements of sigma minus for n = m+3
        Pert_mat_el_offsetbyminusthree = zip_and_add(sigM, Psig)

        # Multiply by perturbative coefficient:
        Pert_mat_el_offsetbyone = list(map(lambda x: Lambda*x, Pert_mat_el_offsetbyone))
        Pert_mat_el_offsetbyminusthree = list(map(lambda x: Lambda*x, \
                                             Pert_mat_el_offsetbyminusthree))

        pert_mat_elems_list.append((mat_el_block.data, \
                                    Pert_mat_el_offsetbyone, \
                                    Pert_mat_el_offsetbyminusthree))

    #    print(list(map(lambda x: list(map(lambda y: shape(y), x)), pert_mat_elems_list)))

    # Now that we have constructed the full set of perturbative matrix
    # elements, we need to construct the perturbative outer products
    # before finally combining them with the matrix elements in order
    # to have the perturbative eigenoperators.

    # Eigenvectors are stored in rows.

    xi_list_nnplusl = []
    for (k, l) in enumerate([1, -1, 3]):
        xi_nnplusl_temp = []
        for i in range(N+1-abs(l)):
            if l > 0:
                ev_list_left, ev_list_right = eig_v.data[i], eig_v.data[i+l]
            else:
                ev_list_left, ev_list_right = eig_v.data[i+abs(l)], eig_v.data[i]
                
            NCi, NCiplusl = len(ev_list_left), len(ev_list_right)
            xi_nnplusl = np.zeros((NCi, NCiplusl, NCi, NCiplusl))

            for m in range(NCi):
                for n in range(NCiplusl):
                    ev_L = ev_list_left[m]
                    ev_R = ev_list_right[n]
                    A = np.outer(ev_L, ev_R)
                    xi_nnplusl[:][:][m][n] = A

            xi_nnplusl_temp.append(xi_nnplusl)

        xi_list_nnplusl.append(xi_nnplusl_temp)

    xi_list_nnplusone = xi_list_nnplusl[0]
    xi_list_nnminusone = xi_list_nnplusl[1]
    xi_list_nnplusthree = xi_list_nnplusl[2]

    # Re-check if following code doesn't work...
    
    # Next, we combine the unperturbed outer products with our
    # perturbative correction matrices in order to form all of the
    # perturbative outer products that we need.

    def zip_and_contract(A, B):
        zipped = zip(A, B)
        contr = map(lambda x: np.einsum('ij, iklm', x[0], x[1]), zipped)
        return list(contr)

    def zip_and_contract_conj(A, B):
        A_conj = map(lambda x: x.conj(), A)
        zipped = zip(A_conj, B)
        contr = map(lambda x: np.einsum('ik, jilm', x[0], x[1]), zipped)
        return list(contr)
    
    Pt_Xi_list = zip_and_contract(pp, xi_list_nnminusone[1:])
    Xi_Mc_list = zip_and_contract_conj(mm, xi_list_nnminusone[:len(xi_list_nnminusone)-1])

    Mt_Xi_list = zip_and_contract(mm[:len(mm)-1], xi_list_nnplusthree)
    Xi_Pc_list = zip_and_contract_conj(pp[1:], xi_list_nnplusthree)

    A_offset_minusone_siteless = [(Pt_Xi_list[0], 0)] + \
                                 list(zip(Pt_Xi_list[1:], Xi_Mc_list[:len(Xi_Mc_list)-1])) + \
                                 [(0, Xi_Mc_list[len(Xi_Mc_list)-1])]
    
    A_offset_plusthree_siteless = [(0, Xi_Pc_list[0]), (0, Xi_Pc_list[1])] + \
                                  list(zip(Mt_Xi_list[:len(Mt_Xi_list)-2], Xi_Pc_list[2:])) + \
                                  [(Mt_Xi_list[len(Mt_Xi_list)-2], 0), (Mt_Xi_list[len(Mt_Xi_list)-1], 0)]

    
    # Now we calculate the full eigenoperators for all the possible
    # spin-changes, weighted by the correct (possibly already
    # perturbed) matric elements
    
    A_l_nnpone = []
    A_l_nnmone = []
    A_l_nnpthree = []
    for pert_mat_elem in pert_mat_elems_list:
        # pert_mat_elem is a triple containing the perturbative
        # matrix element information, for each site.

        mat_el_block = pert_mat_elem[0]
        
        Pert_mat_el_offsetbyone = pert_mat_elem[1]
        Pert_mat_el_offsetbythree = pert_mat_elem[2]

        ###################################################################################
        # n = m+1
        # This is the most complicated case
        A_l_nnplusone_list = []
        for (i, mat_el_offsetbyone) in enumerate(mat_el_block):
            A_l_shape = shape(xi_list_nnplusone[i])
            NCi = A_l_shape[0]
            NCiplusone = A_l_shape[1]
            
            A_l_nnplusone = np.zeros(A_l_shape)

            for j in range(NCi):
                for k in range(NCiplusone):
                    A_l_nnplusone[j][k][:][:] = mat_el_offsetbyone[j][k] * \
                                                 xi_list_nnplusone[i][j][k][:][:]
            A_l_nnplusone_list.append((A_l_nnplusone,
                                       A_offset_minusone_siteless[i],
                                       A_offset_plusthree_siteless[i]))
        A_l_nnpone.append(A_l_nnplusone_list)

        ###################################################################
        # n = m-1
        A_l_nnminusone_list = []
        for (i, mat_el_offsetbyone) in enumerate(Pert_mat_el_offsetbyone):
            A_l_shape = shape(xi_list_nnminusone[i])
            NCi = A_l_shape[0]
            NCiminusone = A_l_shape[1]
            
            A_l_nnminusone = np.zeros(A_l_shape)

            for j in range(NCi):
                for k in range(NCiminusone):
                    A_l_nnminusone[j][k][:][:] = mat_el_offsetbyone[j][k] * \
                                                 xi_list_nnminusone[i][j][k][:][:]
            A_l_nnminusone_list.append(A_l_nnminusone)
        A_l_nnmone.append(A_l_nnminusone_list)

        ####################################################################
        # n = m+3
        A_l_nnplusthree_list = []
        for (i, mat_el_offsetbythree) in enumerate(Pert_mat_el_offsetbythree):
            A_l_shape = shape(xi_list_nnplusthree[i])
            NCi = A_l_shape[0]
            NCiplusthree = A_l_shape[1]
            
            A_l_nnplusthree = np.zeros(A_l_shape)

            for j in range(NCi):
                for k in range(NCiplusthree):
                    A_l_nnplusthree[j][k][:][:] = mat_el_offsetbythree[j][k] * \
                                                 xi_list_nnplusthree[i][j][k][:][:]
            A_l_nnplusthree_list.append(A_l_nnplusthree)
        A_l_nnpthree.append(A_l_nnplusthree_list)            

    # Next we want to put all of the eigenoperators into their
    # respective full 2^N-square matrices in a 1D array, which then
    # may be summed over in all possible ways in order to calculate
    # the complete Liouvillian

    N_C_k = [bin(N, k) for k in range(N+1)]
    N_C_k_Cummulative = [0]
    for k in range(1,N+1):
        N_C_k_Cummulative.append(sum(N_C_k[:k]))
    
    Perturbative_eigenoperators_sites = []
    for l in range(N):
        Perturbative_eigenoperators_l = []
        A_l_nnplusone_tuple = A_l_nnpone[l]
        A_l_nnminusone_list = A_l_nnmone[l]
        A_l_nnplusthree_list = A_l_nnpthree[l]
        
        # For each of these variables, we need to extract the relevant
        # eigenoperators, sum the differently shaped operators in the
        # full 2^N-space and then stick them all in one giant list.

        # Let us first deal with the harder case of
        # A_l_nnplusone_tuple, which has the structure:
        # [(t, ((t, t), (t, t)))]
        
        for (i, t) in enumerate(A_l_nnminusone_list):
            #print(shape(t))
            NCi, NCim1 = shape(t)[0], shape(t)[1]
            for m in range(NCi):
                for n in range(NCim1):
                    omega = eig_en[i+1][m] - eig_en[i][n]
                    Perturbative_eigenoperators_l.append((omega, block_to_2N(t[m][n], N, N_C_k_Cummulative[i], N_C_k_Cummulative[i-1])))
        
        for (i, t) in enumerate(A_l_nnplusthree_list):
            #print(shape(t))
            NCi, NCip3 = shape(t)[0], shape(t)[1]
            for m in range(NCi):
                for n in range(NCip3):
                    omega = eig_en[i][m] - eig_en[i+3][n]
                    Perturbative_eigenoperators_l.append(
                        (omega, block_to_2N(t[m][n], N,
                         N_C_k_Cummulative[i], N_C_k_Cummulative[i+3])))
                                        
        for (i, tp) in enumerate(A_l_nnplusone_tuple):
            e_op = np.zeros((2**N, 2**N))
            current_e_op_list = []
            # let us extract each of the tuples and for each case
            # construct the set of all eigenoperators, before we turn
            # them into full-sized matrices

            # There will be 5 lists of eigenoperators in general,
            # which we need to add together to get the final full
            # perturbed eigenoperator

            sig_xi_nnpone = tp[0]
            NCi, NCiplusone = shape(sig_xi_nnpone)[0], shape(sig_xi_nnpone)[1]
            Pt_xi = tp[1][0]
            Xi_Mc = tp[1][1]
            Mt_xi = tp[2][0]
            Xi_Pc = tp[2][1]
            #print(shape(sig_xi_nnpone), shape(Pt_xi), shape(Xi_Mc), shape(Mt_xi), shape(Xi_Pc))
            
            #print(i)
            for m in range(NCi):
                for n in range(NCiplusone):
                    # now, for each of the five tensors, check if they
                    # exist or not, and then turn them each into full
                    # matrices. Add them at the end and finally append
                    # them to the total eigenoperator count.
                    sparse_mat_list = []
                    sparse_mat_list.append(block_to_2N(sig_xi_nnpone[m][n], N, N_C_k_Cummulative[i], N_C_k_Cummulative[i+1]))
                
                    if not type(Pt_xi) == int:
                        # This one goes backwards from i+2 to i+1
                        NCip2, NCip1 = bin(N, i+2), bin(N, i+1)
                        sparse_mat_list.append(Lambda*block_to_2N(Pt_xi[m][n], N, N_C_k_Cummulative[i+2], N_C_k_Cummulative[i+1]))
                                                
                    if not type(Xi_Mc) == int:
                        # i to i-1
                        NCim1 = bin(N, i-1)
                        sparse_mat_list.append(Lambda*block_to_2N(Xi_Mc[m][n], N, N_C_k_Cummulative[i], N_C_k_Cummulative[i-1]))
                        
                    if not type(Mt_xi) == int:
                        # i-2 to i+1
                        NCim2 = bin(N, i-2)
                        sparse_mat_list.append(Lambda*block_to_2N(Mt_xi[m][n], N, N_C_k_Cummulative[i-2], N_C_k_Cummulative[i+1]))
                        
                    if not type(Xi_Pc) == int:
                        # i to i+3
                        NCip3 = bin(N, i+3)
                        sparse_mat_list.append(Lambda*block_to_2N(Xi_Pc[m][n], N, N_C_k_Cummulative[i], N_C_k_Cummulative[i+3]))

                    #print(shape(sum(sparse_mat_list)))
                    omega = eig_en[i][m] - eig_en[i+1][n]
                    Perturbative_eigenoperators_l.append((omega, sum(sparse_mat_list)))

        Perturbative_eigenoperators_sites.append(Perturbative_eigenoperators_l)
        
    return Perturbative_eigenoperators_sites




# # FUNCTION: pert_master

# DESCRIPTION: Utilise the perturbed eigenoperators to contruct the
#              perturbed Liouvillian.

def master(eops, J, N):

    L = 0.
    ijtot = 0.
    jitot = 0.
    for i, j in eops: 
        ijtot += i.dag()*j
        jitot += j.dag()*i
        L += sprepost(i,j.dag()) + sprepost(j,i.dag())
    L += -spre(ijtot)-spost(jitot)

    return L
