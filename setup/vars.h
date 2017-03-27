integer, parameter :: NFILES = 8
   integer unit_i, unit_o, unit_pot
!
   type(commLineOption) :: files(NFILES) = (/& ! COMPTE: cal que -o sigui el 6 per evitar problemes en debug
      commLineOption("-rlib","reslib","formatted","old","Residue Library"),&
      commLineOption("-pot","potential","formatted","old","DMD Potentials"),&
      commLineOption("-pdbin","pdbin","formatted","old","Input structure(PDB)"),&
      commLineOption("-top","topology","unformatted","unknown","Topology"),&
      commLineOption("-r", "coordinates","unformatted","unknown","Coordinates"), &
      commLineOption("-o","log","formatted","unknown","Log File") , &
      commLineOption("-ligpdbin", "ligpdbin","formatted","old","Input ligand structure(PDB)"), &
      commLineOption("-i", "null","formatted","unknown","Settings") &
      /)
!
   type(residueLibrary) :: resLib
   type(ffprm) :: ff
   type(struc) :: recstr, ligstr, str
   type(stepPotIntDefList) :: stpList
! 
   integer recNatoms
   