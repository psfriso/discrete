!     
! File:   IO.F
! Author: gelpi
!
! Created on 15 de junio de 2012, 13:09
!
!======================================================================
 subroutine readCommLine (files, unit_i, unit_o, NFILES)
 use commLine
!
 integer, intent(IN) :: NFILES
 integer, intent(OUT) :: unit_i, unit_o
 type(commLineOption), intent(INOUT) :: files(NFILES) 
!
 call inputArgs(files)
 unit_i=-1
 if (trim(filename(files, '-i')).ne.'null') unit_i = openFn(files,'-i')
 if (fileName(files,'-o').ne.'log') then
   unit_o = openFn(files, '-o')
 else
   unit_o = 6
 endif
 call header (unit_o)
 end subroutine readCommLine
!======================================================================
 subroutine header (unit_o)
  integer, intent(IN) :: unit_o
   write (unit_o, *) "================================================="
   write (unit_o, *) "=                                               ="
   write (unit_o, *) "=                DMDSetup  (v. 0.2.6)           ="
   write (unit_o, *) "=                   (c) 2012                    ="
   write (unit_o, *) "=                                               ="
   write (unit_o, *) "================================================="
   write (unit_o, *)
 end subroutine header
!======================================================================