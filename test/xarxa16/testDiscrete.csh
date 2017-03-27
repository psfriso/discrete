#!/usr/bin/tcsh
#set number of threads 1..16
set np=1
#
setenv OMP_NUM_THREADS $np
../../discrete/discrete -i dmdtestM.in -sdf xarxa16.sdf -traj xarxa16.$np.traj.pdb -ener xarxa16.$np.ener.txt -rst xarxa16.rst -o $np.out
