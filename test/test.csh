#!/bin/csh
../setup/DMDSet -rlib ../dat/dmd_resLibrary.dat -pot ../dat/dmd_potentials.dat -r $1.r -top $1.top -pdbin $1.pdb 
../discrete/discrete -i dmdtest.in -r $1.r -top $1.top -traj $1.traj.pdb -ener $1.ener.txt -rst $1.rst 
