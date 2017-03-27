#!/bin/csh
foreach f (*pdb)
../../setup/DMDSet -rlib ../../dat/dmd_resLibrary.dat -pot ../../dat/dmd_potentials.dat -r ${f:r}.r -top ${f:r}.top -pdbin $f 
end
