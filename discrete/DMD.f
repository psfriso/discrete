!
! File:   DMD.f
! Author: gelpi
!
! Created on 14 de marzo de 2012, 16:54
!
subroutine dmdIntloop(srt, TACT, TMIN)
    use intList
    use geometryDP
    use Simulation
    use constants
    type(simulRT), intent(INOUT) :: srt
    real*8, intent(IN) :: TACT, TMIN
    !
    real*8 tevent, tevent1, tevent100
    integer mem1, mem2, npair1, npair2, i, j
    !
    srt % tacact = 0.
    srt % toUpdate = .true.
    ! evolucio temporal
    !------------------------------------------------------------------------------
    do while (srt % tacact .lt. TACT)
        ! update colision lists
        do i = 1, srt % natom
            if (srt % toUpdate(i)) then
                srt % ipart(i) = minloc(srt % nblist(i) % iData(1:srt % nblist(i) % nats) % timp, 1)
                srt % tpart(i) = srt % nblist(i) % iData(srt % ipart(i)) % timp
            endif
        enddo
        srt % toUpdate = .false.
        ! busca quina es la propera colisio
        call nextCol(mem1, mem2, npair1, npair2, tevent, srt % ipart, srt % tpart, srt % natom, srt % nblist)
        !
        tevent1 = tevent - srt % temps
        tevent100 = tevent1/100.
        srt % temps = tevent
        srt % tacact = srt % tacact + tevent1
        srt % iev = srt % iev + 1
        ! translacio i variacio dels temps
        do i = 1, srt % natom ! mantenim la versio inline degut a la barreja de tipus de real !!
            srt % r(i) % x = srt % r(i) % x + srt % v(i) % x * tevent1
            srt % r(i) % y = srt % r(i) % y + srt % v(i) % y * tevent1
            srt % r(i) % z = srt % r(i) % z + srt % v(i) % z * tevent1
        enddo
        call updateV(srt % r, srt % v, srt % nblist(mem1) % iData(npair1) % deltak, &
        srt % top % ats(mem1) % xm, srt % top % ats(mem2) % xm, srt % nblist(mem1) % iData(npair1) % xsum, &
        mem1, mem2, srt % natom)
        ! Assegurem que trespassem la barrera en la direcció correcta
        srt % r(mem1) % x = srt % r(mem1) % x + srt % v(mem1) % x * tevent100
        srt % r(mem1) % y = srt % r(mem1) % y + srt % v(mem1) % y * tevent100
        srt % r(mem1) % z = srt % r(mem1) % z + srt % v(mem1) % z * tevent100
        srt % r(mem2) % x = srt % r(mem2) % x + srt % v(mem2) % x * tevent100
        srt % r(mem2) % y = srt % r(mem2) % y + srt % v(mem2) % y * tevent100
        srt % r(mem2) % z = srt % r(mem2) % z + srt % v(mem2) % z * tevent100
        !
        ! ara actualitza els temps de colisio per a les dues particules que han xocat
        ! anulem la colisio per a que no es repeteixi
        srt % nblist(mem1) % iData(npair1) % timp = FSPERS
        srt % nblist(mem2) % iData(npair2) % timp = FSPERS
        call updateXocPart(mem1, srt % nblist, srt % temps, srt % r, srt % v, TMIN, srt % natom, srt % toUpdate)
        call updateXocPart(mem2, srt % nblist, srt % temps, srt % r, srt % v, TMIN, srt % natom, srt % toUpdate)
    enddo
    ! end do while (tacact.lt.tact)------------------------------------------------
end subroutine dmdIntloop
