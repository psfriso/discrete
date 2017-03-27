subroutine colisioBond(srt)
    use geometry
    use geometryDP
    use Simulation
    type(simulRT), intent(INOUT) :: srt
    real csang, calcSignCosAng
    real rbmin, rbmax
    real*8 rij
    integer i, j, k
    do i = 1, srt%natom
        do k = 1, srt%blist(i) % nats
            j = srt%blist(i) % iData(k) % patnum
            csang = calcSignCosAng(srt%r, srt%v, i, j, srt%natom)
            rij = calcDistDP(srt%r(i), srt%r(j))
            rbmin = srt%blist(i) % idata(k) % stepPt % step(1) % r
            rbmax = srt%blist(i) % idata(k) % stepPt % step(2) % r
            if ((rij .gt. rbmax.and.csang .gt. 0.).or.(rij .lt. rbmin.and.csang .lt. 0.)) then
                call chgmom(i, j, srt%top%ats(i)%xm, srt%top%ats(j)%xm, srt%blist(i) % idata(k) % xsum, &
                srt%r, srt%v, srt%natom)
            endif
        enddo
    enddo
    return
end subroutine colisioBond
!===============================================================================
subroutine colisioNonBond(srt, TMIN, ISOLV)
    use geometry
    use geometryDP
    use intList
    use Simulation
    type(simulRT), intent(INOUT) :: srt
    integer, intent(IN) :: ISOLV    
    real, intent(IN) :: TMIN
    real csang, calcSignCosAng
    real*8 rij
    integer i, j, k, l
    !
    ! removes overlaps
    do i = 1, srt%natom
        do k = 1, srt%nblist(i) % nats
            j = srt%nblist(i) % iData(k) % patnum
            if (j .gt. i.and..not.srt%nblist(i) % iData(k) % stepPt % isalt) then
                csang = calcSignCosAng(srt%r, srt%v, i, j, srt%natom)
                rij = calcDist2DP(srt%r(i), srt%r(j))
                if (rij .lt. srt%nblist(i) % iData(k) % stepPt % dmin.and.csang .lt. 0.) then
                    call chgmom(i, j, srt%top%ats(i)%xm, srt%top%ats(j)%xm, srt%nblist(i) % iData(k) % xsum, &
                    srt%r, srt%v, srt%natom)
                    srt%ierr = srt%ierr + 1
                endif
            endif
        enddo
    enddo
    ! calculate colision times
    do i = 1, srt%natom
        do k = 1, srt%nblist(i) % nats
            j = srt%nblist(i) % iData(k) % patnum
            l = srt%nblist(i) % iData(k) % simp
            if (j .gt. i) then
                if (srt%stepPts(i, j) % tipInt .eq. SSEC.or.ISOLV .eq. 1) then
                    call updateTCol(i, srt%temps, srt%r, srt%v, srt%nblist(i) % iData(k), TMIN, srt%natom)
                    srt%nblist(j) % iData(l) % timp = srt%nblist(i) % iData(k) % timp
                    srt%nblist(j) % iData(l) % deltak = srt%nblist(i) % iData(k) % deltak
                endif
            endif
        enddo
    enddo
end subroutine colisioNonBond
!===============================================================================
pure subroutine chgmom(mem1, mem2, xm1, xm2, xsum, r, v, natom)
    use geometry
    use geometryDP
    use Simulation
    real, intent(IN) :: xsum, xm1, xm2
    integer, intent(IN) :: mem1, mem2, natom
    type(pointDP), intent(IN) :: r(natom)

    type(pointDP), intent(INOUT) :: v(natom)

    type(pointDP) dr
    type(pointDP) dv
    real*8 dp
    !
    dr = r(mem2) - r(mem1)
    dv = v(mem2) - v(mem1)
    dp = -(dr % x * dv % x + dr % y * dv % y + dr % z * dv % z) / dotDP(dr, dr) / xsum
    ! modul del moment transferit en la colisio
    v(mem1) % x = v(mem1) % x - dp / xm1 * dr % x
    v(mem1) % y = v(mem1) % y - dp / xm1 * dr % y
    v(mem1) % z = v(mem1) % z - dp / xm1 * dr % z

    v(mem2) % x = v(mem2) % x + dp / xm2 * dr % x
    v(mem2) % y = v(mem2) % y + dp / xm2 * dr % y
    v(mem2) % z = v(mem2) % z + dp / xm2 * dr % z
end subroutine chgmom
!===============================================================================
subroutine updateTCol(i, temps, r, v, intD, TMIN, natom)
    ! input r v de dues particules, output deltak i timp
    use geometry
    use geometryDP
    use stepPotentials
    use intList
    use constants
    integer, intent(IN) :: natom, i
    integer j
    type(pointDP), intent(IN) :: r(natom)
    type(pointDP), intent(IN) :: v(natom)
    type(intData), intent(INOUT) :: intD
    real*8, intent(IN) :: temps
    real, intent(IN) :: TMIN

    integer k
    real*8 argk(MAXSTEPS), tijs(MAXSTEPS, 2)
    real*8 rij, rij2, vij, vij2
    real*8 dotrv, a, b
    type(pointDP) dr
    type(pointDP) dv
    !
    integer lc(2)
    !   
    j = intD % patnum
    intD % timp = FSPERS
    intD % deltak = 0.
    !   
    dr = r(j) - r(i)
    dv = v(j) - v(i)
    rij2 = dotDP(dr, dr)
    rij = sqrt(rij2)
    vij2 = dotDP(dv, dv)
    vij = sqrt(vij2)
    dotrv = dr % x * dv % x + dr % y * dv % y + dr % z * dv % z
    tijs = 1.d15
    if (rij2 .le. intD % stepPt % step(intD % stepPt % nstep) % r**2 + dotrv**2/vij2) then
        do k = 1, intD % stepPt % nstep
            argk(k) = intD % stepPt % step(k) % r**2 - rij2 + dotrv**2/vij2
            if (argk(k) .ge. 0.) then
                a = -dotrv / vij2
                b = sqrt(argk(k)) / vij
                tijs(k, 1) = a - b
                tijs(k, 2) = a + b
            endif
        enddo
        lc = minloc(tijs, mask = tijs .gt. TMIN)
        if (lc(1) .gt. 0.) then
            intD % timp = temps + tijs(lc(1), lc(2))
            if (lc(2) .eq. 1) then
                intD % deltak = intD % stepPt % step(lc(1)) % e
            else
                intD % deltak = -intD % stepPt % step(lc(1)) % e
            endif
        endif
    endif
end subroutine updateTCol
!===============================================================================
subroutine updateV(r, v, deltak, xm1, xm2, xsum, mem1, mem2,natom)
    use geometry
    use geometryDP
    use constants
    integer, intent(IN) :: mem1, mem2,natom
    type(pointDP), intent(IN) :: r(natom)
    real, intent(IN) :: deltak, xsum, xm1, xm2
    type(pointDP), intent(INOUT) :: v(natom)

    type(pointDP) dr
    type(pointDP) dv
    real*8 rij, sto, vdmod, dp
!    real, parameter :: A2 = 1.e-20 * 1.e30

    dr = r(mem2) - r(mem1)
    dv = v(mem2) - v(mem1)
    ! calculo vdmod, la projeccio de la diferencia de velocitats en l'eix que uneix les particules
    rij = dsqrt(dotDP(dr, dr))
    vdmod = (dv % x * dr % x + dv % y * dr % y + dv % z * dr % z) / rij
    ! entra o surt d'un pou
    sto = vdmod**2 + 4. * deltak * xsum / A2
    if (sto .gt. 0) then
        ! vario la velocitat
        dp = -vdmod + sign(1.d0, vdmod) * sqrt(sto)
        dp = dp / 2. / xsum / rij
    else
        ! les particules es queden atrapades al pou
        dp = -vdmod / xsum / rij
    endif
    v(mem1) % x = v(mem1) % x - dp / xm1 * dr % x
    v(mem1) % y = v(mem1) % y - dp / xm1 * dr % y
    v(mem1) % z = v(mem1) % z - dp / xm1 * dr % z

    v(mem2) % x = v(mem2) % x + dp / xm2 * dr % x
    v(mem2) % y = v(mem2) % y + dp / xm2 * dr % y
    v(mem2) % z = v(mem2) % z + dp / xm2 * dr % z
end subroutine updateV
!==============================================================================
subroutine nextCol(mem1, mem2, np1, np2, tevent, ipart, tpart, natom, nblist)
    use intList
    integer, intent(IN) :: natom
    integer, intent(IN) :: ipart(natom)
    real*8, intent(IN) :: tpart(natom)
    type (intpList), intent(IN) :: nblist(natom)
    integer, intent(OUT) :: mem1, mem2
    real*8, intent(OUT) :: tevent
    integer, intent(OUT) :: np1, np2
    !  
    mem1 = minloc(tpart, 1)
    if (mem1 .le. 0.or.ipart(mem1) .le. 0) then
        write (0, *) "ERROR: No colision found (this should never happen!!)"
        stop 1
    endif
    tevent = tpart(mem1)
    np1 = ipart(mem1)
    mem2 = nblist(mem1) % iData(ipart(mem1)) % patnum
    np2 = nblist(mem1) % iData(ipart(mem1)) % simp
end subroutine nextCol
!===============================================================================  
subroutine updateXocPart(m1, nblist, temps, r, v, TMIN, natom, toUpdate)
    use intList
    use geometry
    use geometryDP
    integer, intent(IN) :: natom, m1
    real, intent(IN) :: TMIN
    type(pointDP), intent(IN) :: r(natom)
    type(pointDP), intent(IN) :: v(natom)
    type(intpList), intent(INOUT) :: nblist(natom)
    real*8, intent(IN) :: temps
    integer j, k, l
    logical, intent(INOUT) :: toUpdate(natom)
    !    
    do k = 1, nblist(m1) % nats
        j = nblist(m1) % iData(k) % patnum
        l = nblist(m1) % iData(k) % simp
        call updateTCol(m1, temps, r, v, nblist(m1) % iData(k), TMIN, natom)
        nblist(j) % iData(l) % timp = nblist(m1) % iData(k) % timp
        nblist(j) % iData(l) % deltak = nblist(m1) % iData(k) % deltak
        toUpdate(j) = .true.
    enddo
    toUpdate(m1) = .true.
end subroutine updateXocPart
!===============================================================================   
pure function calcSignCosAng(r, v, i, j, natom)
    use geometryDP
    real calcSignCosAng
    integer, intent(IN) :: natom, i, j
    type(pointDP), intent(IN) :: r(natom)
    type(pointDP), intent(IN) :: v(natom)
    real*8 :: dotv
    type(pointDP) dr
    type(pointDP) dv
    !
    dr = r(j) - r(i)
    dv = v(j) - v(i)
    dotv = dr % x * dv % x + dr % y * dv % y + dr % z * dv % z
    calcSignCosAng = sign(1., real(dotv))
end function calcSignCosAng
!===============================================================================   
