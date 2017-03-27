!===============================================================================   
subroutine readCommLine(files, unit_i, unit_o, NFILES)
    USE commLine
    integer, intent(IN) :: NFILES
    integer, intent(OUT) :: unit_i, unit_o
    type(commLineOption), intent(INOUT) :: files(NFILES)
    !
    call inputArgs(files)
    unit_i = openFn(files, '-i')
    if (fileName(files, '-o') .ne. 'log') then
        unit_o = openFn(files, '-o')
    else
        unit_o = 6
    endif
end subroutine readCommLine
!===============================================================================   
subroutine openFileByName(unitf, fileName, status, fmt)
    integer unitf, err
    character*(*) filename, fmt, status
    open (unitf, file = trim(filename), status = status, form = fmt, iostat = err)
    if (err .ne. 0) then
        write (0, '("Error ",i5, " opening file: ", a50)') err, trim(filename)
        stop 1
    end if
end subroutine openFileByName
!===============================================================================   
subroutine programHeader(unit_o)
    integer unit_o
#ifdef OMP
    integer omp_get_max_threads
#endif
    write (unit_o, *) "================================================="
    write (unit_o, *) "=                                               ="
    write (unit_o, *) "=               DISCRETE  (0.2.7a)              ="
    write (unit_o, *) "=                                               ="
    write (unit_o, *) "=     A. Emperador, J. L. Gelpi, M.Orozco       ="
    write (unit_o, *) "=                                               ="
    write (unit_o, *) "=                  (c) 2012                     ="
    write (unit_o, *) "================================================="
    write (unit_o, *)
#ifdef OMP
    write (unit_o, '(" Running OpenMP version on ", i3, " processors")') omp_get_max_threads()
    write (unit_o, *)
#endif
end subroutine programHeader
!===============================================================================   
subroutine errorAllocmem(ioerr, text)
    integer ioerr
    character(len = *) text
    if (ioerr .ne. 0) then
        write (6, '("Error ",i5, " allocating memory",a30)') ioerr, text
        stop 1
    endif
end subroutine errorAllocmem
!===============================================================================   
subroutine writeSnapshot(unit_traj, ibloc, sysdef, srts, OUTFORMAT)
    use geometryDP
    use Simulation
    integer, intent(IN) :: unit_traj, ibloc, OUTFORMAT
    type(systemDef), intent(IN) :: sysdef
    type(simulRT), intent(in) :: srts(sysdef % nsimul)
    type(atomTopology) at
    integer i, ns

    if (OUTFORMAT .eq. 1) then
        do ns = 1, sysdef % nsimul
            write (unit_traj, '(10f8.3)') (srts(ns) % r(i), i = 1, srts(ns) % natom)
        enddo
    else
        !TODO escriptura repliques
        write (unit_traj, '("MODEL",8X,I4,1X,I4)') ibloc
        do ns = 1, sysdef % nsimul
            write (unit_traj, '("REMARK Simul. No.",8X,I5)') ns
            do i = 1, srts(ns) % natom
                at = srts(ns) % top % ats(i)
                write (unit_traj, '("ATOM",2X,I5,2X,A3,1X,A3,1X,A1,1X,I3,4X,3F8.3)') i, &
                at % atomid, at % res, at % chain, at % rnum, srts(ns) % r(i)
            enddo
            write (unit_traj, '("TER")')
        enddo
    endif
    write (unit_traj, '("ENDMDL")')
end subroutine writeSnapshot
!===============================================================================   
subroutine writeRestartFiles(sysdef, srts, RSTV)
    use geometry
    use geometryDP
    use Simulation
    use txtUtils
    !
    type(systemDef), intent(IN) :: sysdef
    type(simulRT), intent(in) :: srts(sysdef % nsimul)
    integer, intent(IN) :: RSTV
    character(len = 100) nst
    !
    type(point), allocatable :: rsp(:) ! Single precision version of r for binary I/O
    integer unit_rst, unit_rstv, i, ns
    !  
    do ns = 1, sysdef % nsimul
        allocate (rsp(srts(ns) % natom))
        write (nst, '(A50,".",i5,".rst")') sysdef % CoorFiles(srts(ns) % imol), ns
        call removesp(nst)
        call openFileByName(unit_rst, nst, 'unknown', "unformatted")
        do i = 1, srts(ns) % natom
            rsp(i) = DPtoSP(srts(ns)%r(i))
        enddo
        write (unit_rst) srts(ns) % natom
        write (unit_rst) (rsp(i), i = 1, srts(ns) % natom)
        close(unit_rst)
        !
        !    if (RSTV .eq. 1) then
        !        unit_rstv = openFn(files, '-rstv')
        !        write (unit_rstv) natom
        !        write (unit_rstv) (v(i), i = 1, natom)
        !        close(unit_rstv)
        !    endif
        deallocate(rsp)
    enddo
end subroutine writeRestartFiles
!===============================================================================   
subroutine writeEnergies(unit_ener, unit_o, irepl, srt, WENER, WLOG) !pendent parametresps
    use constants
    use Simulation
    type(simulRT), intent(IN) :: srt
    integer, intent(IN) :: unit_ener, unit_o, irepl
    logical, intent(IN) :: WENER, WLOG
    character(len = 20), parameter :: fmtEner = '(1x,i5,10f8.2)'
    character(len = 120), parameter :: fmtLog = &
    '(" R:",i4," TIME: ",f8.2," ps Events: ",i4," Overlaps ",i4," TEMP:" f8.2," K") '
    !    
    character(LEN=*), parameter ::fmtEner2 = '(1X,"EPOT: ",f7.2,1X,7(A,1X,f7.2,1X))'
    integer i
    if (WENER) &
    write (unit_ener, fmtEner) irepl, &
    srt % temps * PSPERFS, &
    srt % epotfis + srt % epotgo, &
    srt % ekin/FACTE, &
    srt % epotfis + srt % epotgo + srt % ekin0/FACTE, &
    srt % ekin0/1.5/srt % natom
    write (unit_ener, fmtEner2) srt%epotfis+srt%epotgo, (etermLab(i),srt % epot(eterm(i)),i=1,NOUTETERMS)
    !
    if (WLOG) then
        write (unit_o, fmtLog) irepl, &
        srt % temps * PSPERFS, &
        srt % iev, &
        srt % ierr, &
        srt % ekin0/1.5/srt % natom
        write (unit_o, fmtEner2) srt%epotfis+srt%epotgo, (etermLab(i), srt % epot(eterm(i)), i = 1, NOUTETERMS)
        write (unit_o, *)
    endif
end subroutine writeEnergies
!===============================================================================   
subroutine writeAllEnergies(unit_ener, unit_o, sysdef, strs, WENER, WLOG)
    use constants
    use Simulation
    type(systemDef), intent(IN) :: sysdef
    type(simulRT), intent(IN) :: strs(sysdef % nsimul)
    integer, intent(IN) :: unit_ener, unit_o
    integer irepl
    logical, intent(IN) :: WENER, WLOG
    do irepl = 1, sysdef % nsimul
        call writeEnergies(unit_ener, unit_o, irepl, strs(irepl), WENER, WLOG)
    enddo
end subroutine writeAllEnergies
!===============================================================================   
