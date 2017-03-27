!
! Subject : Simulation module
! Revision: $id$
!
MODULE Simulation
    use geometry
    use geometryDP
    use intList
    use constants

    type atomTopology
        character(len = 4) atomId, res, atp
        character(len = 1) chain
        integer rnum, oid
        logical dummy, frozen
        real qq, evdw, rvdw, xm, rhc
    end type atomTopology

    type strucTopology
        integer molType, nmol, nres, natom, npots
        type (stepPotIntDef) stepDefList(MAXSPOTDEF)
        type(atomTopology), allocatable :: ats(:)
        type(intDataTop), allocatable :: stepPtsDef(:,:)
        real, allocatable :: xsum(:,:)
        real xmassa
    end type strucTopology

    type simulRef
        integer natom
        type(point), allocatable :: rsp(:) ! Single precision version of r for binary I/O
    end type simulref

    type simulrt
        integer natom
        integer imol ! indice para recuperar la topologia 
        type(stepPotInt), allocatable :: stepPts(:,:)
        type(pointDP), allocatable :: r(:), v(:), f(:)
        type(point) rcm ! C of M   
        real epot(MAXTIPINT), ekin, epotfis, epotgo, ekin0
        real temp
        type(intpList), allocatable :: blist(:), nblist(:)
        logical, allocatable :: toUpdate(:)
        real*8 tacact, taccorr, temps, tacrect
        real*8, allocatable :: tpart(:)
        integer, allocatable :: ipart(:)
        integer iev, ierr
        type(strucTopology), pointer :: top ! pointer to Topology
    end type simulrt

    type systemDef
        integer nmolTypes, nsimul
        integer, allocatable :: numCopies(:)
        character*100, allocatable :: topFiles(:), CoorFiles(:)
        character*50, allocatable :: molId(:)
    end type systemDef
    !===============================================================================
CONTAINS
    !systemDef
    function allocateSystemDef(nmoltypes) result (sysdef)
        type(systemDef) sysdef
        integer, intent(IN) :: nmoltypes
        allocate (sysdef % topFiles(nmolTypes), sysdef % CoorFiles(nmolTypes), &
        sysdef % numCopies(nmolTypes), sysdef % molId(nmolTypes))
        sysdef % nmoltypes = nmoltypes
        sysdef % topFiles = ''
        sysdef % coorFiles = ''
        sysdef % numCopies = 0
        sysdef % molId = ''
    end function allocateSystemDef
    !===============================================================================
    function loadSystemDef(unit) result(sysdef)
        type(systemDef) sysdef
        integer, intent(IN) :: unit
        character*80 str
        character*3 label
        integer i
        i = 0
        100 read(unit, '(a80)', end = 110) str
        if (str(1:3) .eq. 'MOL') i = i + 1
        goto 100
        110 rewind(unit)
        sysdef = allocateSystemDef(i)
        sysdef % nmoltypes = i
        i = 0
        120 read(unit, '(a80)', end = 130) str
        if (str(1:3) .eq. 'MOL') then
            i = i + 1
            read (str, '(a3,i5,a50)') label, sysdef % numCopies(i), sysDef % molId(i)
            read (unit, *) sysdef % topFiles(i)
            read (unit, *) sysdef % coorFiles(i)
        endif
        goto 120
        130 continue
        sysdef % nsimul = sum(sysdef % numCopies)
    end function loadSystemDef
    !===============================================================================   
    !strucTopology
    subroutine allocateStrucTopology(top)
        type (strucTopology), intent(INOUT) :: top
        integer ioerr
        allocate (top % ats(top % natom), &
        top % xsum(top % natom, top % natom), &
        top % stepPtsDef(top % natom, top % natom), stat = ioerr)
        if (ioerr .ne. 0) call errorAllocmem(ioerr, 'StructTopology')
    end subroutine allocateStrucTopology
    !===============================================================================
    function loadTopology(unit_top) result (top)
        use constants
        type(strucTopology) top
        integer, intent(IN) :: unit_top
        integer itv, i, j
        !    
        read(unit_top) itv
        if (itv .ne. itopVersion) then
            write(0, *) "Error: incompatible topology, ", itopVersion, " required, got ", itv
            stop 1
        endif
        read(unit_top) top % molType
        read(unit_top) top % nmol, top % nres, top % natom
        call allocateStrucTopology(top)
        !
        read (unit_top) (top % ats(i) % atomId, top % ats(i) % rnum, top % ats(i) % res, &
        top % ats(i) % chain, top % ats(i) % atp, top % ats(i) % frozen, top % ats(i) % dummy, i = 1, top % natom)
        !
        do i = 1, top % natom
            read(unit_top) top % ats(i) % qq, top % ats(i) % evdw, top % ats(i) % rvdw, top % ats(i) % rhc, top % ats(i) % xm
        enddo
        !
        read (unit_top) top % npots
        !      write (unit_o, '( i4, " Step Potentials definitions loaded")') top%npots
        do i = 1, top % npots
            read (unit_top) top % stepDefList(i)
        enddo
        !    
        100 read (unit_top, end = 110) i, j, top % stepPtsDef(i, j)
        goto 100
        110 continue
        do j = 2, top % natom
            do i = 1, j - 1
                top % xsum(i, j) = 0.5 * (1./top % ats(i) % xm + 1./top % ats(j) % xm)
                top % xsum(j, i) = top % xsum(i, j)
            enddo
        enddo
    end function loadTopology
    !===============================================================================   
    !simulRef
    function allocateSimulRef(natom) result (srt)
        integer, intent(IN) :: natom
        type(simulref) :: srt
        integer ioerr
        srt % natom = natom
        allocate (srt % rsp(natom), stat = ioerr)
        if (ioerr .gt. 0) call errorAllocmem(ioerr, 'allocateSimulRef')
    end function allocateSimulRef
    !===============================================================================
    function loadCoordinates(unit_r, natom) result(refStruc)
        type(simulRef) refStruc
        integer, intent(IN) :: unit_r, natom
        integer i, j
        read (unit_r) j
        if (j .ne. natom) then
            write (0, *) "ERROR: coordinates and topology files do not match"
            stop 1
        endif
        refStruc = allocateSimulRef(natom)
        read (unit_r) (refStruc % rsp(i) % x, refStruc % rsp(i) % y, refStruc % rsp(i) % z, i = 1, natom)
    end function loadCoordinates
    !===============================================================================   
    !simulRT
    function allocateSimulrt(natom) result (srt)
        integer, intent(IN) :: natom
        type(simulrt) :: srt
        integer i, ioerr
        srt % natom = natom
        allocate ( &
        srt % stepPts(natom, natom), &
        srt % r(natom), &
        srt % v(natom), &
        srt % f(natom), &
        srt % blist(natom), &
        srt % nblist(natom), &
        srt % tpart(natom), &
        srt % ipart(natom), &
        srt % toUpdate(natom), stat = ioerr)
        if (ioerr .gt. 0) call errorAllocmem(ioerr, 'allocateSimulRT')
        do i = 1, natom
            srt % blist(i) = allocateintPList(natom, ioerr)
            if (ioerr .gt. 0) call errorAllocmem(ioerr, 'allocateSimulRT - BList')
            srt % nblist(i) = allocateintPList(natom, ioerr)
            if (ioerr .gt. 0) call errorAllocmem(ioerr, 'allocateSimulRT - NBList')
        enddo
    end function allocateSimulrt
    !===============================================================================
    function cloneSimulRT(str) result(newstr)
        type(simulRT), intent(IN) :: str
        type(simulRT) newstr
        newstr = allocateSimulRT(str % natom)
        newStr = str
    end function cloneSimulRT
    !===============================================================================
    function minExt(Srt) result (p)
        use geometryDP
        type(simulRT), intent(IN) :: srt
        type(pointDP) p
        p % x = minval(srt % r % x)
        p % y = minval(srt % r % y)
        p % z = minval(srt % r % z)
    end function minExt
    !===============================================================================
    function maxExt(Srt) result (p)
        use geometryDP
        type(simulRT), intent(IN) :: srt
        type(pointDP) p
        p % x = maxval(srt % r % x)
        p % y = maxval(srt % r % y)
        p % z = maxval(srt % r % z)
    end function maxExt
    !===============================================================================
    function center(Srt) result (p) ! Geometric
        use geometryDP
        type(simulRT), intent(IN) :: srt
        type(pointDP) p, p1, p2
        p1 = minExt(Srt)
        p2 = maxExt(Srt)
        p = p1 + p2
        p = prodFactDP(0.5d0, p)
    end function center
    !===============================================================================
    subroutine calcCM(srt) ! center of mass
        use geometry
        use geometryDP
        type(simulRT), intent(INOUT) :: srt
        type(point) rcm
        integer i
        !
        srt % rcm = point(0., 0., 0.)
        do i = 1, srt % natom
            srt % rcm % x = rcm % x + srt % top % ats(i) % xm * real(srt % r(i) % x)
            srt % rcm % y = rcm % y + srt % top % ats(i) % xm * real(srt % r(i) % y)
            srt % rcm % z = rcm % z + srt % top % ats(i) % xm * real(srt % r(i) % z)
        enddo
        srt % rcm = (1./srt % top % xmassa) * srt % rcm
    end subroutine calcCM
    !===============================================================================
    subroutine moveToOrigin(srt)
        use geometryDP
        type(simulRT), intent(INOUT) :: srt
        integer i
        do i = 1, srt % natom
            srt % r(i) = srt % r(i) - SPtoDP(srt % rcm)
        enddo
    end subroutine moveToOrigin
    !===============================================================================   
    subroutine prepareSimStruc(newTop, reftop, newStr, Str)
        type (strucTopology), intent(IN) :: reftop
        type (strucTopology), intent(INOUT) :: newTop
        type (simulRef), intent(IN) :: str
        type (simulRT), intent(INOUT) :: newStr
        integer i, j, nat
        newTop % natom = refTop % natom - count(refTop % ats % dummy)
        call allocateStrucTopology(newTop)
        newTop % stepDefList = refTop % stepDefList
        newStr = allocateSimulRT(newTop % natom)
        nat = 0
        do i = 1, refTop % natom
            if (.not.reftop % ats(i) % dummy) then
                nat = nat + 1
                newTop % ats(nat) = reftop % ats(i)
                newTop % ats(nat) % oid = i
                newStr % r(nat) = sPtoDP(str % rsp(i))
            endif
        enddo
        do j = 2, newTop % natom
            do i = 1, j - 1
                newTop % stepPtsDef(i, j) = refTop % stepPtsDef(newTop % ats(i) % oid, newTop % ats(j) % oid)
                newTop % xsum(i, j) = refTop % xsum(newTop % ats(i) % oid, newTop % ats(j) % oid)
                newTop % xsum(j, i) = newTop % xsum(i, j)
            enddo
        enddo
        newTop % xmassa = sum(newTop % ats % xm)
    end subroutine prepareSimStruc

END MODULE Simulation
!===============================================================================
