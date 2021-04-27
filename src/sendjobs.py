import subprocess
import sys
import os

import numpy as np

# python sendjob.py Nmin Nmax Nstep gmin gmax gstep omega0min omega0max omega0step bmin bmax bstep
# Where omega0 is in units of bandwidth, which is precalculated

[Nmin, Nmax, Nstep] = [int(sys.argv[i]) for i in range(1,4)]
[gmin, gmax, gstep] = [float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])]
[omega0min, omega0max, omega0step] = [float(sys.argv[7]), float(sys.argv[8]), int(sys.argv[9])]
[bmin, bmax, bstep] = [float(sys.argv[10]), float(sys.argv[11]), int(sys.argv[12])]

print([Nmin, Nmax, Nstep])
print([gmin, gmax, gstep])
print([omega0min, omega0max, omega0step])
print([bmin, bmax, bstep])


def if_exist_move(directory):
    if os.path.isdir(directory):
        os.chdir(directory)
    else:
        os.makedirs(directory)
        os.chdir(directory)
    return 0

for N in range(Nmin, Nmax, Nstep):
    Ndir = 'N='+str(N)
    if_exist_move(Ndir)

    for g in np.linspace(gmin, gmax, gstep):
        gdir = 'g='+str(g)
        if_exist_move(gdir)

        for omega in np.linspace(omega0min, omega0max, omega0step):
            omdir = 'omega0='+str(omega)+'xBandWidth'
            if_exist_move(omdir)

            for b in np.linspace(bmin, bmax, bstep):
                bdir = 'b='+str(b)
                if_exist_move(bdir)
                fname = 'script.sh'

                with open(fname, 'w') as script:
                    script.write('#!/bin/bash\n')
                    script.write('#$ -S /bin/bash\n')
                    script.write('#$ -cwd\n')
                    script.write('#$ -pe smp 1\n')
                    script.write('#$ -j y\n')
                    script.write('export LD_LIBRARY_PATH=~/miniconda3/lib\n')
                    script.write('~/miniconda3/bin/python ~/Calculations/xymodel/xy_structured_bath.py '\
                                     +str(N)+' '+str(g)+' '+str(omega)+' '+str(b))

                subprocess.call(['qsub', fname])
                os.chdir('..')

            os.chdir('..')

        os.chdir('..')

    os.chdir('..')

