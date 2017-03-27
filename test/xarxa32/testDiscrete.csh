#!/usr/bin/tcsh
#set number of threads 1..32
set np=1
#
setenv OMP_NUM_THREADS $np
../../discrete/discrete -i dmdtestM.in -sdf xarxa16.sdf -traj xarxa16.$np.traj.pdb -ener xarxa32.$np.ener.txt -rst xarxa32.rst -o $np.out
