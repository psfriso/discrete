subroutine buildStepPotentials(srt, IDAB, rcutcoul2, rcutsolv2)
    use simulation
    use constants
    use geometryDP
    type(simulRT), intent(INOUT) :: srt
    integer IDAB
    real rcutcoul2, rcutsolv2
    integer i, j
    real distat2(srt % natom, srt % natom)
    !
    do j = 2, srt % natom
        do i = 1, j - 1
            distat2(i, j) = real(calcDist2DP(srt % r(i), srt % r(j)))
            distat2(j, i) = distat2(i, j)
        enddo
    enddo
    do j = 2, srt % natom
        do i = 1, j - 1
            if (srt % top % stepPtsDef(i, j) % id .gt. 0) then
                srt % stepPts(i, j) = getStepPotFromDEF(&
                srt % top % stepDefList(srt % top % stepPtsDef(i, j) % id), &
                srt % top % stepPtsDef(i, j) % rref, &
                srt % top % stepPtsDef(i, j) % eref, .false., &
                srt % top % stepPtsDef(i, j) % tipInt)
                srt % stepPts(i, j) % dmin = (srt % top % ats(i) % rhc + srt % top % ats(j) % rhc)**2
                srt % stepPts(i, j) % isalt = (srt % top % ats(j) % rnum .eq. srt % top % ats(i + 1) % rnum.and.(&
                (srt % top % ats(i) % atomid .eq. 'C' .and. srt % top % ats(j) % atomid .eq. 'CB').or. &
                (srt % top % ats(i) % atomid .eq. 'CB'.and. srt % top % ats(j) % atomid .eq. 'N') .or. &
                (srt % top % ats(i) % atomid .eq. 'O' .and. srt % top % ats(j) % atomid .eq. 'CA')))
            endif
        enddo
    enddo
    do j = 2, srt % natom
        do i = 1, j - 1
            !SSEC
            if (srt % stepPts(i, j) % tipint .eq. SSEC) then
                if (distat2(i, j) .gt. RCUTGO2.or.IDAB.eq.0) then
                    srt % stepPts(i, j) % tipInt = NULL
                    srt % stepPts(i, j) % active = .false.
                else
                    srt % stepPts(i, j) % active = .true.
                endif
            ! BOND terms
            elseif (srt % stepPts(i, j) % tipInt .eq. COVB.or. &
                srt % stepPts(i, j) % tipInt .eq. COVF.or. &
                srt % stepPts(i, j) % tipInt .eq. HB) then
                srt % stepPts(i, j) % active = .false.
                ! j > i only 
                srt % blist(i) % nats = srt % blist(i) % nats + 1
                srt % blist(i) % iData(srt % blist(i) % nats) = intData(j, srt % blist(j) % nats, &
                srt % stepPts(i, j), srt % top % xsum(i, j), FSPERS, 0.)
            endif
        enddo
    enddo
    call activateStepPot(srt, rcutcoul2, rcutsolv2)
end subroutine buildStepPotentials
!===============================================================================
subroutine activateStepPot(srt, rcutcoul2, rcutsolv2)
    use Simulation
    use geometryDP
    use constants
    type(simulRT), intent(INOUT) :: srt
    real, intent(IN) :: rcutcoul2, rcutsolv2
    integer i, j
    real*8 rij2
    !
    !SSEC is kept active. Pending: Variable SSEC restrains
    where (srt % stepPts % active.and.srt%stepPts%tipInt.ne.SSEC)
        srt % stepPts % active = .false.
    end where
    srt % nblist % nats = 0
    do j = 2, srt % natom
        do i = 1, j - 1
            if (srt % stepPts(i, j) % tipInt .ne. SSEC) then
                rij2 = calcDist2DP(srt % r(i), srt % r(j))
                srt % stepPts(i, j) % active = &
                ((rij2 .lt. rcutcoul2.and.srt % stepPts(i, j) % tipInt .eq. COUL).or. &
                (rij2 .lt. rcutsolv2.and.srt % stepPts(i, j) % tipInt .eq. HFIL).or. &
                (rij2 .lt. rcutsolv2.and.srt % stepPts(i, j) % tipInt .eq. HFOB)) 
            endif
        enddo
    enddo
    do j = 2, srt % natom
        do i = 1, j - 1
            if (srt%stepPts(i,j)%active) then
                srt % nblist(j) % nats = srt % nblist(j) % nats + 1
                srt % nblist(i) % nats = srt % nblist(i) % nats + 1
                srt % nblist(i) % iData(srt % nblist(i) % nats) = intData(j, srt % nblist(j) % nats, &
                srt % stepPts(i, j), srt % top % xsum(i, j), FSPERS, 0.)
                srt % nblist(j) % iData(srt % nblist(j) % nats) = intData(i, srt % nblist(i) % nats, &
                srt % stepPts(i, j), srt % top % xsum(i, j), FSPERS, 0.)
            endif
        enddo
    enddo
end subroutine activateStepPot
!========================================================================
subroutine calcEpot(srt)
    use geometryDP
    use stepPotentials
    use constants
    use Simulation

    type(simulRT), intent(INOUT) :: srt
    real dist
    integer i, j, k
    !PENDENT TREBALLAR AMB NBLIST
    srt % epotgo = 0.
    srt % epotfis = 0.
    srt % epot = 0.
    do j = 2, srt % natom
        do i = 1, j - 1
            if (srt % stepPts(i, j) % active) then
                dist = sqrt(real(calcDist2DP(srt % r(i), srt % r(j))))
                k = srt % stepPts(i, j) % nstep
                do while ((k .gt. 1).and.dist .lt. srt % stepPts(i, j) % step(k) % r)
                    srt % epot(srt % stepPts(i, j) % tipInt) = srt % epot(srt % stepPts(i, j) % tipInt) - &
                    srt % stepPts(i, j) % step(k) % e / FACTE
                    k = k - 1
                enddo
                if (dist .lt. srt % stepPts(i, j) % step(k) % r) &
                srt % epot(srt % stepPts(i, j) % tipInt) = srt % epot(srt % stepPts(i, j) % tipInt) - &
                srt % stepPts(i, j) % step(k) % e / FACTE
            endif
        enddo
    enddo
    srt % epotgo = srt % epot(SSEC)
    srt % epotfis = sum(srt % epot) - srt % epotgo
end subroutine calcEpot
!========================================================================
pure function calcEkin(srt) result (ekin)
use geometryDP
use constants
use Simulation
type(simulRT), intent(IN) :: srt
real ekin
integer i
ekin = 0.
do i = 1, srt % natom
    ekin = ekin + 0.5 * srt % top % ats(i) % xm * real(A2 * dotDP(srt % v(i), srt % v(i)))
enddo
end function calcEkin
!===============================================================================
subroutine thermalize(seed, srt, RMVCM, TCALC, TCORR)
    use geometryDP
    use Simulation
    type(simulRT), intent(INOUT) :: srt
    integer, intent(IN) :: TCALC, TCORR, RMVCM
    integer, intent(IN) :: seed
    if (TCALC .gt. 0) then
        if (TCALC .eq. 1) then
            call rescaleV(seed, srt, RMVCM)
        else if (TCALC .eq. 2) then
            call andersenBath(seed, srt, TCORR, 0) !RMVCM = 0 mentre no hi pensem
        endif
    endif
end subroutine thermalize
!===============================================================================
subroutine rescaleV(seed, srt, RMVCM)
    use geometryDP
    use Simulation
    use random
    type(simulRT), intent(INOUT) :: srt
    integer, intent(IN) :: seed, RMVCM
    integer i, kk
    real calcEkin
    type(pointDP) vcm
    real*8 fi, sto
    !
    vcm = pointDP(0., 0., 0.)
    do i = 1, srt % natom
        kk = -(seed + 3 * i + 1 + srt % iev)
        fi = ran1(kk)
        srt % v(i) % x = fi - 0.5
        kk = -(seed + 3 * i + 2 + srt % iev)
        fi = ran1(kk)
        srt % v(i) % y = fi - 0.5
        kk = -(seed + 3 * i + 3 + srt % iev)
        fi = ran1(kk)
        srt % v(i) % z = fi - 0.5
        !
        vcm = vcm + srt % top % ats(i) % xm * 1.d0 * srt % v(i)
    enddo
    if (RMVCM .eq. 1) then
        vcm = (1.d0/srt % top % xmassa) * vcm
        do i = 1, srt % natom
            srt % v(i) = srt % v(i) - vcm
        enddo
    endif
    sto = sqrt(1.5 * srt % natom * srt % TEMP / calcEkin(srt))
    do i = 1, srt % natom
        srt % v(i) = sto * srt % v(i)
    enddo
    srt % ekin = 1.5 * srt % natom * srt % TEMP ! No cal recalcularla 
end subroutine rescaleV
!===============================================================================
subroutine andersenBath(seed, srt, TCORR, RMVCM)
    use geometryDP
    use Simulation
    use random
    !
    type(simulRT), intent(INOUT) :: srt
    integer, intent(IN) :: TCORR, RMVCM
    integer, intent(IN) :: seed
    integer i, kk
    real calcEkin
    type(pointDP) vcm
    real*8 fi, sto, ff
    !
    ! Andersen thermostat
    !
    vcm = pointDP(0., 0., 0.)
    kk = -(seed + srt % iev)
    fi = ran1(kk)
    i = int(srt % natom * fi) + 1
    ff = sqrt(srt % TEMP/srt % top % ats(i) % xm)/1.d5
    kk = -(seed + 3 * i + 1 + srt % iev)
    fi = ran1(kk)
    srt % v(i) % x = (2. * fi - 1.) * ff
    kk = -(seed + 3 * i + 2 + srt % iev)
    fi = ran1(kk)
    srt % v(i) % y = (2. * fi - 1.) * ff
    kk = -(seed + 3 * i + 3 + srt % iev)
    fi = ran1(kk)
    srt % v(i) % z = (2. * fi - 1.) * ff
    !
    do i = 1, srt % natom
        srt % v(i) = srt % v(i) + TCORR/srt % top % ats(i) % xm * 1.d0 * srt % f(i)
    enddo
!    if (RMVCM .eq. 1) then
!        vcm = (1.d0/srt % top % xmassa) * vcm
!        do i = 1, srt % natom
!            srt % v(i) = srt % v(i) - vcm
!        enddo
!    endif
    sto = sqrt(1.5 * srt % natom * srt % TEMP / calcEkin(srt))
    do i = 1, srt % natom
        srt % v(i) = sto * srt % v(i)
    enddo
    srt % ekin = 1.5 * srt % natom * srt % TEMP ! No cal recalcularla 
end subroutine andersenBath
!========================================================================
subroutine updateForces(sysdef, srts)
    use Simulation
    type(systemDef), intent(IN) :: sysdef
    type(simulRT) srts(sysdef % nsimul)
    integer i, j
    do i = 1, sysdef % nsimul - 1
        do j = i + 1, sysdef % nsimul
            call forceij(srts(i) % top, srts(j) % top, srts(i), srts(j))
        enddo
    enddo
end subroutine updateForces
!========================================================================
subroutine forceij(topi, topj, srti, srtj)
    use geometryDP
    use Simulation
    type(strucTopology), intent(IN) :: topi, topj
    type(simulRT), intent(INOUT) :: srti, srtj
    type(pointDP) dr
    real*8 rij, rij2, fcoul, fvdw
    real qij, epsij, sigmaij
    integer i, j
    do i = 1, topi % natom
        do j = 1, topj % natom
            ! forces entre particules de diferents subsistemes
            dr = srtj % r(j) - srti % r(i)
            rij2 = calcdist2DP(srti % r(i), srtj % r(j))
            rij = sqrt(rij2)
            qij = topi % ats(i) % qq * topj % ats(j) % qq
            fcoul = qij/rij**3
            ! coulomb
            srti % f(i) % x = srti % f(i) % x + fcoul * dr % x
            srti % f(i) % y = srti % f(i) % y + fcoul * dr % y
            srti % f(i) % z = srti % f(i) % z + fcoul * dr % z
            !
            srtj % f(j) % x = srtj % f(j) % x - fcoul * dr % x
            srtj % f(j) % y = srtj % f(j) % y - fcoul * dr % x
            srtj % f(j) % z = srtj % f(j) % z - fcoul * dr % x
            !! Van der Waals
            epsij = sqrt(topi % ats(i) % evdw * topj % ats(j) % evdw)
            sigmaij = topi % ats(i) % rvdw + topj % ats(j) % rvdw
            fvdw = epsij * (sigmaij**12/rij**14 - 2. * sigmaij**6/rij**8)
            srti % f(i) % x = srti % f(i) % x + fvdw * dr % x
            srti % f(i) % y = srti % f(i) % y + fvdw * dr % y
            srti % f(i) % x = srti % f(i) % z + fvdw * dr % z
            !
            srtj % f(j) % x = srtj % f(j) % x - fvdw * dr % x
            srtj % f(j) % y = srtj % f(j) % y - fvdw * dr % y
            srtj % f(j) % x = srtj % f(j) % z - fvdw * dr % z
        enddo
    enddo
end subroutine forceij
!========================================================================
