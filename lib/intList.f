 MODULE intList
 use stepPotentials
 
 type intData
    integer patnum, simp ! patnum: num atom de la parella, simp: index en llista de la parella
    type(stepPotInt) :: stepPt
    real :: xsum
    real*8 :: timp
    real :: deltak    
 end type intData
 
 type intpList
    integer nats
    type(intData), pointer :: iData(:)
 end type intpList
 
 CONTAINS
 
 function allocateintPList(natom, ioerr) result (pl)
 integer, intent(IN) :: natom
 integer, intent(OUT) :: ioerr
 type (intpList) pl
  allocate (pl%iData(natom), stat=ioerr)
  pl%nats=0
 end function allocateintPList
 
 function cloneIntList(natom,lst) result(newlst)
 type(intpList), intent(IN) :: lst
 type(intPlist) newlst
 integer, intent(IN) :: natom
 integer ioerr,i
 newlst = allocateintPList(natom, ioerr)
 newLst%nats=lst%nats
 do i=1,newLst%nats
     newLst%iData(i)=lst%iData(i)
 enddo
 end function cloneIntList
 END MODULE intList
