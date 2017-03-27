!
! Authors : Josep Ll. Gelpi, Agusti Emperador
! Subject : Main Setup program for DMD
! Revision: $id$
!
! v. 0.2. Prepared for multiple molecules separated by TER, possible coverage of EDMD
! v, 0.2.1. Prepared for protein-protein docking. implements multiscale en frozen atoms
! v. 0.2.6. Adds setup for potentials, preparation for CG, new topology format. Version number to match discrete
!
program DMDSetup
use Constants
use paramSet
use commLine
use Structure
use potentials
use resLibrary
use stepPotentials
!
#include "vars.h"
   integer i,k,l
   !
  call readCommLine (files, unit_i, unit_o, NFILES)
  call printFileSummary(files, unit_o)
   write (unit_o,*)
!
   if (unit_i.gt.0) then
      call readInputParamSet(unit_i)
      call printInputParamSet(unit_o)
   endif
!      
  resLib = loadResLibrary(openFn(files, '-rlib'))
   write(unit_o,'(1X, "Residue Library: " i4 " residues / ", i5, " atoms loaded")') resLib%nres, sum(resLib%residues%natoms)
  unit_pot = openFn(files, '-pot')
  ff = loadPotentials(unit_pot)
  write(unit_o,'(1X, "Atom Types:      ", i4, " types")') ff%ntypes
   write(unit_o,*)
  stpList = loadStPotIntDef(unit_pot)
   write(unit_o,'(1X, "StepPotentials: ", i4," types")') stpList%npots
   write (unit_o,*)
   
   if (TIPCALC.eq.MD) then
     str = loadStructurePDB(openFn(files, '-pdbin'))
      write(unit_o,'(1X, "Structure:       ", i4, " molecules / ", i4, " residues / ", i5, " atoms loaded")') &
            str%nmol, str%nres, str%natoms
      write(unit_o, '(" Molecule ",i3,": ",2i5)') (i, str%mols(i),i=1,str%nmol)
      recNatoms=str%natoms
  elseif (TIPCALC.eq.DOCKING) then
     write(unit_o,'(1X,"Setup for Protein-protein docking")') 
      recstr = loadStructurePDB(openFn(files, '-pdbin'))
      write(unit_o,'(1X, "Receptor Structure:    ", i4, " molecules / ", i4, " residues / ", i5, " atoms loaded")') &
      recstr%nmol, recstr%nres, recstr%natoms
      write(unit_o, '(" Receptor Molecule",i3,": ",2i5)') (i, recstr%mols(i),i=1,recstr%nmol)
      recNatoms=recStr%natoms
!
    ligstr = loadStructurePDB(openFn(files, '-ligpdbin'))
      write(unit_o,'(1X, "Ligand Structure:      ", i4, " molecules / ", i4, " residues / ", i5, " atoms loaded")') &
         ligstr%nmol, ligstr%nres, ligstr%natoms
      write(unit_o, '(" Ligand Molecule  ",i3,": ",2i5)') (i, ligstr%mols(i),i=1,ligstr%nmol)
      write(unit_o,'(" Building complex...")')         
     str = mergeStructures (recstr, ligstr, offsetx, offsety,offsetx)         
     write(unit_o,'(1X, "Complex Structure:     ", i4, " molecules / ", i4, " residues / ", i5, " atoms loaded")') &
         str%nmol, str%nres, str%natoms
      write(unit_o, '(" Complex Molecule  ",i3,": ",2i5)') (i, str%mols(i),i=1,str%nmol)
      write(unit_o,*)
!
     call markInterface(str, recNatoms, dint0, dint1, irig, SetCoreCa, SetFrozen)
      write (unit_o,'(1X,"Getting interface...")')
      write (unit_o,'(1X,"Interface dist. ", f5.1," All Atom layer: ",f5.1)') dint0, dint1
      write (unit_o,'(1X,"           Receptor    Ligand")')
      write (unit_o,'(1X,"Interface  ",         i7,"    ", i7)') count(.not.str%ats(1:recStr%natoms)%frozen), &
         count(.not.str%ats(recStr%natoms+1:str%natoms)%frozen)
      write (unit_o,'(1X,"Frozen     ",         i7,"    ", i7)') & 
         count(str%ats(1:recStr%natoms)%frozen)-count(str%ats(1:recStr%natoms)%dummy), &
         count(str%ats(recStr%natoms+1:str%natoms)%frozen)-count(str%ats(recStr%natoms+1:str%natoms)%dummy) 
      write (unit_o,'(1X,"Discarded  ",         i7,"    ", i7)') count(str%ats(1:recStr%natoms)%dummy), &
         count(str%ats(recStr%natoms+1:str%natoms)%dummy)
      write (unit_o, *)         
   else
      write(0,'(" Calculation type not implemented")')
      stop 1
  endif
            
  call assignAtType(str,resLib)
  call setBonds(str, resLib, tipCalc, recNatoms, bondCutoff)

   write (unit_o,*) 'Hydrogen bonds' 
   do i=1,str%hbs%npairs
      k = str%hbs%list(i,1)
      l = str%hbs%list(i,2)
      write (unit_o,'(2X,i4," ",a15,": ",i4," ",a15,f8.3)')  & 
         str%ats(k)%molnum, writeAtomId(str%ats(k)), & 
         str%ats(l)%molnum, writeAtomId(str%ats(l)), &
         sqrt(str%distat2(k,l))
   enddo
!   
  call calcSecStr(str)
   
   write (unit_o, *) str%nhelix, " Alfa Helices"
   do i=1, str%nhelix
      write (unit_o,'(1x,2i3," ",a10," - ",i3, " ", a10)') i, &
         str%res(str%helix(i)%ini)%molres, writeResidueId(str%res(str%helix(i)%ini)), &
         str%res(str%helix(i)%fin)%molres, writeResidueId(str%res(str%helix(i)%fin))
   enddo
 !
   write (unit_o, *) str%nbeta, " Beta strands"
   do i=1, str%nbeta
      write (unit_o,'(1x,2i3," ",a10," - ",i3, " ",a10)') i,  &
         str%res(str%beta(i)%ini)%molres, writeResidueId(str%res(str%beta(i)%ini)), & 
         str%res(str%beta(i)%fin)%molres, writeResidueId(str%res(str%beta(i)%fin))
   enddo
!
  call setStepPotentials(str,stpList,ff,fvdw,fsolv)
!   
  call saveTopology(openFn(files,'-top'),str,stpList,ff)
  call saveCoords(openFn(files,'-r'),str)
   write (unit_o,*) "Topology & Coodinates files saved"
 end
 
