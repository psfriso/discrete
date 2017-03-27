module random
        integer, parameter :: b8 = selected_real_kind(14)
contains
    function ran1(idum)
        implicit none !note after use statement
        real(b8) ran1
        integer, intent(inout), optional :: idum
        real(b8) r(97), rm1, rm2
        integer, parameter :: m1 = 259200, ia1 = 7141, ic1 = 54773
        integer, parameter :: m2 = 134456, ia2 = 8121, ic2 = 28411
        integer, parameter :: m3 = 243000, ia3 = 4561, ic3 = 51349
        integer j
        integer iff, ix1, ix2, ix3
        data iff /0/
        save ! corrects a bug in the original routine

        if (present(idum))then
            if (idum < 0.or.iff .eq. 0)then
                rm1 = 1.0_b8/m1
                rm2 = 1.0_b8/m2
                iff = 1
                ix1 = mod(ic1 - idum, m1)
                ix1 = mod(ia1 * ix1 + ic1, m1)
                ix2 = mod(ix1, m2)
                ix1 = mod(ia1 * ix1 + ic1, m1)
                ix3 = mod(ix1, m3)
                do j = 1, 97
                    ix1 = mod(ia1 * ix1 + ic1, m1)
                    ix2 = mod(ia2 * ix2 + ic2, m2)
                    r(j) = (real(ix1, b8) + real(ix2, b8) * rm2) * rm1
                enddo
                idum = 1
            endif
        endif
        ix1 = mod(ia1 * ix1 + ic1, m1)
        ix2 = mod(ia2 * ix2 + ic2, m2)
        ix3 = mod(ia3 * ix3 + ic3, m3)
        j = 1 + (97 * ix3)/m3
        if (j > 97.or.j < 1)then
            write(*, *) ' error in ran1 j=', j
            stop
        endif
        ran1 = r(j)
        r(j) = (real(ix1, b8) + real(ix2, b8) * rm2) * rm1
        return
    end function ran1
END MODULE random




!!  
!!  rutina ran1
!!
!   subroutine ran1(z, idum)
!   integer iff, m1, m2, m3, ia1, ia2, ia3
!   real*8 rm1, rm2, r
!   real*8 z
!   integer ix1, ix2, ix3, j, ic1, ic2, ic3, idum
!!
!   SAVE
!!
!   dimension r(97)
!   parameter (m1 = 259200, ia1 = 7141, ic1 = 54773, rm1 = 1. / m1)
!   parameter (m2 = 134456, ia2 = 8121, ic2 = 28411, rm2 = 1. / m2)
!   parameter (m3 = 243000, ia3 = 4561, ic3 = 51349)
!!
!   data iff /0/
!!
!   if(idum.lt.0 .or. iff.eq.0) then
!      iff = 1
!      ix1 = mod(ic1 - idum , m1)
!      ix1 = mod(ia1 * ix1 + ic1, m1)
!      ix2 = mod(ix1, m2)
!      ix1 = mod(ia1 * ix1 + ic1, m1)
!      ix3 = mod(ix1, m3)
!      do j = 1, 97
!         ix1 = mod(ia1 * ix1 + ic1, m1)
!         ix2 = mod(ia2 * ix2 + ic2, m2)
!         r(j) = (float(ix1) + float(ix2) * rm2) * rm1
!      end do
!   end if
!!
!!   except when initializing, this is were we start.
!!
!   ix1 = mod ( ia1 * ix1 + ic1 , m1 )
!   ix2 = mod ( ia2 * ix2 + ic2 , m2 )
!   ix3 = mod ( ia3 * ix3 + ic3 , m3 )
!   j = 1 + (97 * ix3) / m3
!!
!   if (j.gt.97.or.j.lt.1) write(0,*) ' AAAUUGGHH !!!'
!!
!   z = r(j)
!   r(j) = (float(ix1) + float(ix2) * rm2) * rm1
!!
!   return
!   end
