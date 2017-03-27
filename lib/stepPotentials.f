MODULE stepPotentials
    use Constants
!    integer, parameter :: MAXSTEPS = 3, MAXSPOTDEF = 20
!    real, parameter :: VDWBAR = 10.
!    integer, parameter :: &
!        WELL = 1, & !Well centered in a position. 3 Param DELTA, and EMIN, TODO: infWall .true means always rebound ! 
!        FIT = 2, & !Fits curve A/r^n Param A, n, rcut (last step), nsteps. Step 0 es always rvdwij, use eref to scale
!        USER = 3 ! User provided: Nsteps, r and e arrays (1..nsteps) use eref to scale, 0 step at rvdwij
!    integer, parameter :: NULL=0, COVB = 1, SSEC = 2, HB=3, COUL = 4, HFIL = 5, HFOB = 6, COVF = 7

    type stepPot
        real r, e
    end type stepPot

    type stepPotInt
        integer nstep, tipInt
        type (stepPot) step(MAXSTEPS)
        logical active, isalt
        real dmin
    end type stepPotInt

    type stepPotIntDef
        character(len = 4) :: id
        integer nstep, stepTip
        real r(MAXSTEPS), e(MAXSTEPS)
    end type stepPotIntDef
    
    type intDataTop
        integer :: id, tipInt
        real :: rref, eref
    end type intDataTop    
    
    type stepPotIntDefList
        integer npots
        type(stepPotIntDef) list(MAXSPOTDEF)
    end type stepPotIntDefList

CONTAINS
    !===============================================================================
    function getStepPotFromDEF(stDef, rref, eref, active, tipInt) result (st)
        type(stepPotIntDef), intent(IN) :: stDef
        type(stepPotInt) :: st
        real, intent(IN) :: rref, eref
        logical, intent(IN) :: active
        integer, intent(IN) :: tipInt
        integer i
        st%nstep=0
        select case (stDef % stepTip)
        case (WELL) ! r(1) = DELTA, e(1) = EGO
            st % nstep = 2
            st % step(1) = stepPot((1. - stDef % r(1)) * rref, -stDef % e(1) * FACTE)
            st % step(2) = stepPot((1. + stDef % r(1)) * rref, stDef % e(1) * FACTE)
        case (FIT) ! A/r^n  r(1) = rmax, e(1) = n
            st % nstep = stDef % nstep + 1
            st % step(1) = stepPot(rref, VDWBAR*FACTE*abs(eref))
            do i = 1, int((stDef % r(1) - rref)/stDef % nstep)
                st % step(i) % r = i * (stDef % r(1) - rref)/stDef % nstep + rref
            enddo
            do i = 1, int((stDef % r(1) - rref)/stDef % nstep)
                st % step(i) % e = eref /(st % step(i) % r + st % step(i + 1) % r/2)**stDef % e(1) * FACTE
            enddo
        !TODO Increments en %e
        case (USER)
            st % nstep = stDef % nstep
!            st % step(1) = stepPot(rref, VDWBAR*eref*FACTE)
            do i = 1, st % nstep
                st % step(i) = stepPot(stDef % r(i), stDef % e(i) * eref * FACTE)
            enddo
           st%step(1)%r = rref
        case default
        end select
        st % tipInt = tipInt
        st % active = active
        st % isalt = .false.
        st % dmin = 0.
!        call writeStPotDef(6,stDef)
!        call writeStPot(6,st)
    end function getStepPotFromDef

  !===============================================================================
    function getStepSSec(SIGMAGO, EGO, dist, active) result (st)
        type (stepPotInt) st
        type (stepPotIntDef) stDef
        real SIGMAGO, EGO, dist
        logical active
        stDef = stepPotIntDef ( 'SS',1,WELL,(/SIGMAGO,0.,0./),(/EGO,0.,0./))
        st = getStepPotFromDEF(stDef, dist, 0., active, SSEC)
    end function getStepSSec
!===============================================================================
    function getStepCoul(rvdwij, DPSINT, DPSEXT, ecoul, HPS, active) result (st)
        type( stepPotInt) st
        type (stepPotIntDef) stDef
        real rvdwij, DPSINT, DPSEXT, ecoul, HPS
        logical active
        stDef = stepPotIntDef ( 'COUL',3,USER,(/0.,DPSINT,DPSEXT/),(/-sign(3., ecoul), HPS-1.,-HPS,0./));
        st = getStepPotFromDEF(stDef, rvdwij, ecoul, active, COUL)
    end function getStepCoul
!===============================================================================
    function getStepCoulAuto(rvdwij, DPSEXT, ecoul,active) result (st)
        type( stepPotInt) st
        type (stepPotIntDef) stDef
        real rvdwij,DPSEXT, ecoul
        logical active
        stDef = stepPotIntDef ( 'COULA',3,FIT,(/DPSEXT,0.,0./),(/1., 0.,0./));
        st = getStepPotFromDEF(stDef, rvdwij, ecoul, active, COUL)
    end function getStepCoulAuto
!===============================================================================
    function getStepSolv(rvdwij, DHF, esolv, active) result (st)
        type( stepPotInt) st
        type (stepPotIntDef) stDef
        real rvdwij, DHF, esolv
        logical active
        if (esolv.lt.0.) then
            stDef = stepPotIntDef ('HFOB',2,USER,(/0., DHF,0./),(/3.,-1.,0./));
        else
            stDef = stepPotIntDef ('HFIL',2,USER,(/0.,DHF,0./),(/-1.5,-1.,0./));
        endif
        st = getStepPotFromDEF(stDef, rvdwij, esolv, active, HFOB)
    end function getStepSolv
!===============================================================================
    subroutine writeStPotDef (unit,stDef)
        type(stepPotIntDef) stDef
        integer unit, i
        write (unit,'(1X,A4,2X,2i5)') stDef%id, stDef%stepTip, stDef%nstep
        do i=1,stDef%nstep
            write (unit,'(2F10.4)') stDef%r(i), stDef%e(i)
        enddo
    end subroutine writeStPotDef
!===============================================================================
    subroutine writeStPot (unit,st)
        type(stepPotInt) st
        integer i, unit
        write (unit,'(1X,2I5,2L,F7.3)') st%tipInt, st%nstep, st%active, st%isalt, st%dmin
        do i=1,st%nstep
            write (unit,'(2F10.4)') st%step(i)%r, st%step(i)%e/FACTE
        enddo
    end subroutine writeStPot
!===============================================================================
!TODO
!    subroutine drawStPot(unit,st)
!        type(stepPotInt) st
!        integer unit
!        real iniE,minE,maxE,maxR
!        write (unit,'(1X,2I5,2L,F7.3)') st%tipInt, st%nstep, st%active, st%isalt, st%dmin
      !  iniE = sum(st%step(1:st%nstep)%e)/FACTE
      !  minE = minVal(st%step(1:st%nstep)%e)/FACTE
      !  maxE = maxVal(st%step(1:st%nstep)%e)/FACTE
      !  maxE = maxE + 5.
      !  maxR = maxVal(st%step(1:st%nstep)%r)
      !  do i=0,maxR,0.5
      !  enddo
!    end subroutine drawStPot
 !===============================================================================          
    function loadStPotIntDef(unit) result (stList)
        type(stepPotIntDefList) stList
        type(stepPotIntDef) stDef
        integer unit, i
        character(len = 80) str
        !        
        stList % npots = 0
        10 read(unit, '(a80)', end = 20) str
               write (6, *) str
        if (str(1:4) .eq. 'ENDS') goto 20
        if (str(1:1) .ne. ' '.and.str(1:1) .ne. '#') then
            read (str, '(a4,2i2,A40)', err = 20) stDef % id, stDef % stepTip, stDef % nstep
            do i = 1, stDef % nstep
                read (UNIT, *) stDef % r(i), stDef % e(i)
            enddo
            stList % npots = stList % npots + 1
            stList % list(stList % npots) = stDef
        endif
        goto 10
        20 continue
    end function loadStPotIntDef
!===============================================================================          
    function getstepPotDef (stL, id) result(stDef)
         type(stepPotIntDefList) stL
        type(stepPotIntDef) stDef
        character(*) id
        integer i
        i=1
        do while (trim(stL%list(i)%id).ne.trim(id).and.i.lt.stL%npots)
            i=i+1
        enddo
        if (trim(stL%list(i)%id).eq.trim(id)) then
            stDef = stL%list(i)
        else
            write(0,*) "Requested Potential ", id, " not found"
            stop 1
        endif
    end function getstepPotDef
!===============================================================================        
    function getstepPotDefIndex (stL, id) result(i)
        type(stepPotIntDefList) stL
        character(*) id
        integer i
        i=1
        do while (trim(stL%list(i)%id).ne.trim(id).and.i.lt.stL%npots)
            i=i+1
        enddo
        if (trim(stL%list(i)%id).ne.trim(id)) then
            write(0,*) "Requested Potential ", id, " not found"
            stop 1
        endif
    end function getstepPotDefIndex
!===============================================================================        
    END MODULE stepPotentials
