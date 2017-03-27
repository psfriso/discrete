 MODULE paramSet
     use stepPotentials
     use constants
! Input param
!
   real, save :: &
      temp = 300., &
      rcutgo = 10., &
      rcutcoul = 10., &
      rcutsolv = 7., &
      tmin = 1.e-22*FSPERS, &
      urnomax=RNOMAX, &
      urnomin=RNOMIN, &
      urcomax=RCOMAX, &
      urcomin=RCOMIN, &
      rcutcoul2=0., &
      rcutsolv2=0., &
      urcutgo2=0., &
      xbeta=1000., &
      sclim=0.1
  integer, save :: &
      isolv = 1, &
      nbloc = 1000, &
      idab = 1, &
      igoab = 1, &
      iwr = 1, &
      seed = 2381, &
      TCALC = 1, & ! E Const 0,  T Const 1, TBany Andersen 2
      rmVCM = 1, & ! Force V CM = 0
      outformat = 0, &
      idims = 0, &
      rst=0, &
      rstv=0, &
      nrepl = 1, &
      tiprepl = 0, & ! None, 1 Crowding, 2: Replica Exchange, ...
      multipletop = 0, &
      keepCoords = 1, & ! keep original coordinates, 0:generate random positions for CM's TODO
      replica = 0, & ! Replicated simulations TODO
      splittraj = 0, &
      writelog = 0 ! Log key=value
  real*8, save :: &
      tsnap = 1000., &
      tcorr = 100., &
      tact = 50., &
      trect = 100., &
      tini = 0.
CONTAINS
!===============================================================================
 subroutine readInputParamSet (unit_i)
   integer, intent(IN) :: unit_i
!   
   namelist /input/ tsnap,tcorr,tact,temp,seed,&
            nbloc,rcutgo,tmin,&
            isolv,idab,igoab,iwr,rcutsolv,rcutcoul,&  
            urnomax,urnomin,urcomax,urcomin,rcutgo, tcalc, outformat, trect, xbeta, sclim, idims, &
            rst, tini, rstv, nrepl, multipletop, writelog
!
   read (unit_i, INPUT)
   IDIMS=0 ! forced, DIMS moved to a different executable
   NREPL = max(NREPL,1)
   ! checking 
   if (TRECT.gt.TSNAP) TRECT = TSNAP
   if (TCORR.gt.TRECT) TCORR = TRECT
   if (TACT.gt.TCORR)  TACT = TCORR
   urcutgo2 = rcutgo**2
   rcutcoul2 = rcutcoul**2
   rcutsolv2 = rcutsolv**2
   if (RST.eq.1) RSTV = 1
 end subroutine readInputParamSet
!===============================================================================
 subroutine writeInputParamSet (unit_o)
   integer, intent(IN) :: unit_o
   ! Copiem l'arxiu de parametres. Pendent format
   write (unit_o, *)
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | CALCULATION PARAMETERS                                   |")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Simulation settings                                      |")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Simulation Time (ps) (Nbloc x TSnap)       |",f12.3," |")') NBLOC * TSNAP / 1.e3
   write (unit_o, '(" | Output structure (fs)             | TSnap  |",f12.3," |")') TSNAP 
   if (IDIMS.eq.1) &
   write (unit_o, '(" | Re-scoring target (fs)            | Trect  |",f12.3," |")') TRECT
   write (unit_o, '(" | Update velocities (fs)            | Tcorr  |",f12.3," |")') TCORR
   write (unit_o, '(" | Update Lists, collision times (fs)| Tact   |",f12.3," |")') TACT
   write (unit_o, '(" | Min. accepted colision time (fs)  | TMin   |",f12.8," |")') TMIN   
   write (unit_o, '(" | Temperature (K)                   | Temp   |",6X,f6.2, " |")') TEMP
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | CutOff Beta restrains (A)         | RCutGo |",7X,f5.2, " |")') RCUTGO
   write (unit_o, '(" | Electrostatic cutoff (A)          |RcutCoul|",7X,f5.2, " |")') sqrt(rcutcoul2)
   write (unit_o, '(" | Solvation cutoff (A)              |RcutSolv|",7X,f5.2, " |")') sqrt(rcutsolv2)
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Other                                                    |")')  
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Random generator seed                      |",7X,i5  " |")') seed
   write (unit_o, '(" | IDAB, IGOAB, IWR, ISOLV                    |",4X,4i2," |")') IDAB, IGOAB, IWR, ISOLV
   write (unit_o, '(" ------------------------------------------------------------")')
 end subroutine writeInputParamSet
!===============================================================================
 END MODULE paramSet