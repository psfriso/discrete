!===========================================================================
!  DISCRETE 0.2.7a
!
!  Discrete Molecular Dynamics,
!
!   v. 0.1a    First public release.
!      Read Topology and coordinates prepared by auxiliary program
!        Restart available
!   v. 0.2a    Fortran Free Form.
!      Dynamic Memory access
!   v. 0.2.3a  Optimization memory access, implemented ISOLV=0
!   v. 0.2.4a  Multilevel representation + DIMS
!   v. 0.2.4.1a Multilevel representation (DIMS moved to new code). Vels as A/fs, time as fs
!               Energy conservation bug fixed
!   v. 0.2.6a   Transferencia setup potencials a setup, modificacio format top, preparacio per CG
!   v. 0.2.7a   New data structure for simulation, allow multiple parallel simulations
!   v. 0.2.8a   Test MPI version
!==========================================================================
program discrete
    use commLine
    use geometry
    use geometryDP
    use stepPotentials
    use intList
    use paramSet
    use Simulation
    !
#include "vars.h"
#ifdef OMP
    integer omp_get_max_threads
#endif
    !
    call cpu_time(tinit)
    call readCommLine(files, unit_i, unit_o, NFILES)
    call programHeader(unit_o)
    call printFileSummary(files, unit_o)
    call readInputParamSet(unit_i)
    !
    call writeInputParamSet(unit_o)
    !
    if (MULTIPLETOP .eq. 1) then
        sysDef = loadSystemDef(openFn(files, '-sdf'))
        write (unit_o, *) "System Definition loaded"
        write( unit_o, '(1X,"Num. Mol. Types: ",i5)') sysdef % nmoltypes
        do i = 1, sysdef % nmolTypes
            write (unit_o, '(1X,"Mol. Num: ",i5," ",a)') i, trim(sysDef % molId(i))
            write (unit_o, '(1X,"  -Topology:    ", a)') trim(sysDef % topFiles(i))
            write (unit_o, '(1X,"  -Coordinates: ", a)') trim(sysDef % coorFiles(i))
            write (unit_o, '(1X,"  -Num. copies: ", i5)') sysDef % numCopies(i)
        enddo
    else
        sysDef = allocateSystemDef(1)
        sysDef % numCopies(1) = NREPL
        sysDef % topFiles(1) = fileName(files, '-top')
        sysDef % CoorFiles(1) = fileName(files, '-r')
        sysDef % molId = ''
        sysDef % nsimul = NREPL
    endif
    write (unit_o, '(1X,"Num. simulations:   ", i5)') sysDef % nsimul
    allocate(refTop(sysdef % nmolTypes), refStruc(sysdef % nmolTypes))
    write (unit_o, *)
    write (unit_o, '(" Reading system from topology file(s) ")')
    write (unit_o, *)
    do i = 1, sysDef % nmolTypes
        write (unit_o, *) "Loading ", trim(sysdef % topFiles(i))
        call openFileByName(unit_top, sysdef % topFiles(i), 'old', "unformatted")
        refTop(i) = loadTopology(unit_top)
        write (unit_o, '(" Molecules: "i4," Residues:   ",i4," Atoms:  ",i6)') refTop(i) % nmol, refTop(i) % nres, refTop(i) % natom
        write (unit_o, *) "Loading ", trim(sysdef % coorFiles(i))
        call openFileByName(unit_r, sysdef % CoorFiles(i), 'old', "unformatted")
        refStruc(i) = loadCoordinates(unit_r, refTop(i) % natom)
        close(unit_top)
        close(unit_r)
        write (unit_o, *)
    enddo
    write (unit_o, '(" Reference topologies and coordinates read in ")')
    write (unit_o, *)
    ! Preparing simulations
    allocate(simTop(sysDef % nmoltypes), simStruc(sysDef % nsimul))
    write (unit_o, *) " Allocating ", sysDef%nsimul, " simulations for ", sysdef%nmolTypes, " mol. types"
    write (unit_o,*)
    write (unit_o, '(" Initial Pairs list")')
    write (unit_o, '(" ==================")')
    
    do i = 1, sysDef % nmolTypes
        call prepareSimStruc(simTop(i), reftop(i), simStruc(i), refStruc(i))
        simStruc(i) % top => simTop(i)
        simStruc(i) % imol = i
        call buildStepPotentials(simStruc(i), IDAB, rcutcoul2, rcutsolv2)
        !        
        write (unit_o, '(" Molecule:           ",i5," ",a)') i, trim(sysdef % molId(i))
        write (unit_o, '(" Total:              ",i9)') simTop(i) % natom * (simTop(i) % natom - 1) / 2
        write (unit_o, '(" Bonded:             ",i9)') count(simStruc(i) % stepPts % tipInt .eq. COVB)+ &
        count(simStruc(i) % stepPts % tipint .eq. COVF) + count(simStruc(i) % stepPts % tipInt .eq. HB)
        write (unit_o, '(" Non Bonded:         ",i9)') count(simStruc(i) % stepPts % active)
        write (unit_o, '("   Secondary Struc.  ",i9)') count(simStruc(i) % stepPts % tipInt .eq. SSEC)
        write (unit_o, '("   Electrostatic     ",i9)') count(simStruc(i) % stepPts % tipInt .eq. COUL.and. &
        simStruc(i) % stepPts % active)
        write (unit_o, '("   Hydrophilic       ",i9)') count(simStruc(i) % stepPts % tipInt .eq. HFil.and. &
        simStruc(i) % stepPts % active)
        write (unit_o, '("   Hydrophobic       ",i9)') count(simStruc(i) % stepPts % tipInt .eq. HFob.and. &
        simStruc(i) % stepPts % active)
        !
        call calcEpot(simStruc(i))
        call calcCM(simStruc(i))
        !moveToORigin(simStruc(i))
    enddo
    !
    write (unit_o, *)
    write (unit_o, '(" System setup completed")')
    write (unit_o, *)
    ! TOFO Deallocate
    simStruc % iev = 0
    simStruc % temp = TEMP
    simStruc % temps = 0.
    ns = sysDef % nMolTypes
    do i = 1, sysdef % nMolTypes
        do j = 2, sysdef % numCopies(i)
            ns = ns + 1
            simStruc(ns) = cloneSimulRT(simStruc(i))
        enddo
    enddo
    write (unit_o, '(i5," parallel simulations initialized")') ns
    write (unit_o,*) 
    !    
    if (REPLICA .eq. 0) then
        if (KEEPCOORDS .eq. 0) then
            !TODO Inicializacion coordenadas de moleculas replicadas
        endif
    else
        !TODO inicialización detalles de las replicas TEMP?
    endif
    do i = 1, sysDef % nsimul
        call thermalize(seed, simStruc(i), RMVCM, 1, TCORR)
        simStruc(i) % ekin0 = simStruc(i) % ekin
    enddo
    !    
    !Writing initial energies and coordinates
    !
    write (unit_o, '(" Initial energy evaluation completed")')
    unit_ener = openFn(files, '-ener')
    do i = 1, sysDef % nsimul
        write (unit_o, '(" R:", i4, " EPOT: ",f7.2," EKIN: ",f7.2)') &
        i, simStruc(i) % epotfis, simStruc(i) % ekin/FACTE
        write (unit_ener, '("#R:",i4," Initial energy ",4f10.3)') &
        i, simStruc(i) % epotfis, simStruc(i) % epotgo, simStruc(i) % ekin/FACTE, &
        simStruc(i) % epotfis + simStruc(i) % epotgo + simStruc(i) % ekin/FACTE
        write (unit_ener, '("# time, epot, ekin, Etot, Temp")')
    enddo
    call writeAllEnergies(unit_ener, unit_o, sysdef, simStruc, .true., .false.)
    if (OUTFORMAT .eq. 1) then
            unit_pdb = openFn(files, '-traj')
            unit_traj = openFn(files, '-x')
        call writeSnapshot(unit_pdb, ibloc, sysdef, simStruc, 0)
            write (unit_traj, '(a132)') "(10f8.3)"
        else
        unit_traj = openFn(files, '-traj')
        call writeSnapshot(unit_traj, ibloc, sysdef, simStruc, OUTFORMAT)
        endif
    call cpu_time(tsetup)
#ifdef OMP
    write (unit_o, '("Running test OMP on ",i3," processor(s)")') omp_get_max_threads()
#endif
    write (unit_o,*)
    write (unit_o, '("Energies in kcal/mol")')
    !----------------------------------------------------------------------------
    do ibloc = 1, NBLOC
        tacum = 0.
        do while (tacum .lt. TSNAP)
            simStruc % tacrect = 0.
            !$omp parallel do private(irepl)
            do irepl = 1, sysDef % nsimul
                !TODO inici simulacio replica
                !----------------------------------------------------------------------------
                do while (simStruc(irepl) % tacrect .lt. TRECT)
                    simStruc(irepl) % taccorr = 0.
                    simStruc(irepl) % iev = 0
                    simStruc(irepl) % ierr = 0
                    !----------------------------------------------------------------------------
                    do while (simStruc(irepl) % taccorr .lt. TCORR)
                        call colisioBond    (simStruc(irepl))
                        call activateStepPot(simStruc(irepl), rcutcoul2, rcutsolv2)
                        call colisioNonBond (simStruc(irepl), TMIN, ISOLV)
                        call DMDIntLoop     (simStruc(irepl), TACT, TMIN)
                        simStruc(irepl) % taccorr = simStruc(irepl) % taccorr + simStruc(irepl) % tacact
                    enddo
                    ! end do while(taccorr.lt.tcorr)-----------------------------------------------
                    simStruc(irepl) % ekin0 = calcEkin(simStruc(irepl))
                    call thermalize(seed, simStruc(irepl), RMVCM, TCALC, TCORR)
                    simStruc(irepl) % tacrect = simStruc(irepl) % tacrect + simStruc(irepl) % taccorr
                enddo
                ! end DO while(tacrect.lt.trect) -------------------------------------------------
                if (TCALC .gt. 1) call updateForces(sysDef, simStruc)
                !TODO END repliques
            enddo
            tacum = tacum + maxval(simStruc % tacrect)
        enddo
        ! energies
        do irepl = 1, sysDef % nsimul
            call calcEpot(simStruc(irepl))
        enddo
        call writeAllEnergies(unit_ener, unit_o, sysdef, simStruc, .true., .true.)
        ! trajectory snapshot        
        call writeSnapshot(unit_traj, ibloc, sysdef, simStruc, OUTFORMAT)
        ! restart
        call writeRestartFiles(sysdef, simStruc, RSTV)
        ! end DO while(tacum.lt.tsnap) -------------------------------------------------
    enddo
    ! end DO iblock-----------------------------------------------------------------
    close(5)
    call cpu_time(tfin)
    write (unit_o, *)
    write (unit_o, '(" T I M I N G S ")')
    write (unit_o, '(" ============= ")')
    write (unit_o, '(" Setup: " f10.2," s")') tsetup - tinit
    write (unit_o, '(" Traj:  ",f10.2," s (",f10.2," ns/h, ",f7.2," ns/day)")') tfin - tsetup, &
    simStruc(1) % temps / (tfin - tsetup) * 3.6 / 1000., simStruc(1) % temps / (tfin - tsetup) * 3.6 / 1000. * 24
    write (unit_o, *)
end
!=============================================================================