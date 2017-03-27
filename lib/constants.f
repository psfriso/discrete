!     
! File:   constants.f
! Author: gelpi
!
! Created on 20 de julio de 2012, 8:30
!

MODULE constants
    !sizes
    integer, parameter :: MAXATYPES = 50 ! Atom types
    integer, parameter :: MAXRES = 50, MAXATPERRES = 30 ! residue libraries
    integer, parameter :: MAXSTEPS = 3, MAXSPOTDEF = 20 ! Step Potentials

    !StepPotentials
    integer, parameter :: &
    WELL = 1, & !Well centered in a position. 3 Param DELTA, and EMIN, TODO: infWall .true means always rebound ! 
    FIT = 2, & !Fits curve A/r^n Param A, n, rcut (last step), nsteps. Step 0 es always rvdwij, use eref to scale
    USER = 3 ! User provided: Nsteps, r and e arrays (1..nsteps) use eref to scale, 0 at step 0 means rvdwij
    integer, parameter :: &
    NULL = 0, & ! No interaction
    COVB = 1, & ! Actual Covalent bond
    SSEC = 2, & ! Secondary structure restrain
    HB = 3, & ! Hydrogen Bond
    COUL = 4, & ! Electrostatic
    HFIL = 5, & ! Solvation Hydrophilic
    HFOB = 6, & ! Solvation Hydrophobic
    COVF = 7 ! Covalent restrain forced
    integer, parameter :: MAXTIPINT = 7
    !Output energy
    integer, parameter :: NOUTETERMS = 4
    integer, parameter :: eterm(NOUTETERMS) = (/SSEC, COUL, HFIL, HFOB/)
    character(len=*), parameter :: etermLab(NOUTETERMS) = (/'SSEC: ', 'COUL: ', 'HFIL: ', 'HFOB: '/)
    ! energy
    real, parameter :: FACTE = 4167. ! Energy conversion
    real, parameter :: VDWBAR = 10. ! VDW Wall for automatic step potentials
    !!distances
    real, parameter :: RSSMAX = 2.5
    real, parameter :: RNOMAX = 4.1
    real, parameter :: RNOMIN = 2.5
    real, parameter :: RNCMAX = 5.
    real, parameter :: RNCMIN = 3.2
    real, parameter :: RCOMAX = 5.
    real, parameter :: RCOMIN = 3.2
    real, parameter :: RCUTNB2 = 30. * 30.
    real, parameter :: RCUTGO2 = 20. * 20.
    !structure
    integer, parameter :: PROT = 1, NUC = 2, SMALL = 3, COMPLEX = 4 ! molType
    integer, parameter :: ALL = 0, HEAVY = 1, CAONLY = 2 ! CGTYpe 
    integer, parameter :: HELIX = 1, BETA = 2

    integer, parameter :: MD = 0, DOCKING = 1 ! tipCalc
    integer, parameter :: itopVersion = 0261
    integer, parameter :: PDBFORMAT=0, CRDFORMAT=1
    !
    character(len = 50), parameter :: pdbinputFmt = '(13X,A3,1X,A4,A1,I4,4X,3F8.3,2f6.2)'
    real, parameter :: MINREAL = 1.e-20
    real*8, parameter :: A2 = 1.d-20*1.d30
    real*8, parameter :: FSPERS=1.e15
    real, parameter :: PSPERFS=1.e-3

CONTAINS

END MODULE constants
