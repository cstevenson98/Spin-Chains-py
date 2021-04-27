# A script to make all of the stagerred integrated susceptibility plots for the 6 data points:

import getopt
import sys
import os
from qutip import *
from matplotlib import *

try:
    opts, args = getopt.getopt(sys.argv[1:], "a:p:u", ["PBC", "upd"])
except getopt.GetoptError:
    print('python make_ISS_plots.py -a N,g,omega0,b -p(u)(upd) lam [--PBC]')


flag = 'Unperturbed'
filehead = "rhoSS-bruteforce-N="
boundary_conditions = 'OBC'
Lambda = 0.

for opt, arg in opts:
    if opt == '-h':
        print('python make_ISS_plots.py -a <N,g,omega0,b>')
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
        filehead = "rhoSS::PERT::N="
    elif opt in ("--upd"):
        flag = 'UnperturbedDissipator'
        Lambda = float(arg)
        filehead = "rhoSS::U-PDiss::N="
    elif opt == "--PBC":
        boundary_conditions = 'PBC'

# Construct the sigmax operators for increasing i
Sx_list = []
for i in range(N):
    idd_N = [identity(2) for j in range(N)]
    idd_N[i] = sigmay()
    Sx_list.append(tensor(idd_N))

    
for I in range(6):
    if I + 1 == 1:
        #
        my_list = [N, abs_g, abs_omega0, b, flag, Lambda]
        my_list = list(map(str, my_list))
        filename = filehead+my_list[0]+"-g="+my_list[1]+"-omega0="\
                   +my_list[2]+"-b="+my_list[3]
        
    elif I + 1 == 2:
        #
        my_list = [N, abs_g, 0.0, b, flag, Lambda]
        my_list = list(map(str, my_list))
        filename = filehead+my_list[0]+"-g="+my_list[1]+"-omega0="\
                   +my_list[2]+"-b="+my_list[3]
        
    elif I + 1 == 3:
        #
        my_list = [N, abs_g, -abs_omega0, b, flag, Lambda]
        my_list = list(map(str, my_list))
        filename = filehead+my_list[0]+"-g="+my_list[1]+"-omega0="\
                   +my_list[2]+"-b="+my_list[3]
        
    elif I + 1 == 4:
        #
        my_list = [N, -abs_g, abs_omega0, b, flag, Lambda]
        my_list = list(map(str, my_list))
        filename = filehead+my_list[0]+"-g="+my_list[1]+"-omega0="\
                   +my_list[2]+"-b="+my_list[3]
        
    elif I + 1 == 5:
        #
        my_list = [N, -abs_g, 0.0, b, flag, Lambda]
        my_list = list(map(str, my_list))
        filename = filehead+my_list[0]+"-g="+my_list[1]+"-omega0="\
                   +my_list[2]+"-b="+my_list[3]
        
    elif I + 1 == 6:
        #
        my_list = [N, -abs_g, -abs_omega0, b, flag, Lambda]
        my_list = list(map(str, my_list))
        filename = filehead+my_list[0]+"-g="+my_list[1]+"-omega0="\
                   +my_list[2]+"-b="+my_list[3]
        

    rhoSS = qload(filename)
    
    # For each of the loaded density matrices, calculate the staggered
    # integrated susceptibility as a function of position for
    # increasing i, so S_0,S_i.

    print(I)
    stag_int_susc = 0.
    for k in range(N-1):
        stag_int_susc += ((-1)**(k)) * expect(Sx_list[k]*Sx_list[k+1], rhoSS)
        print(stag_int_susc)
    stag_int_susc += ((-1)**(N-1)) * expect(Sx_list[N-1]*Sx_list[0], rhoSS)
    print(stag_int_susc)
    print('------------')
