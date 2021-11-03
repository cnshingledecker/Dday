PROGRAM test
  CHARACTER(LEN=10), ALLOCATABLE         :: array1(:)
  CHARACTER(LEN=10), ALLOCATABLE         :: array2(:)
  
  ALLOCATE(array1(3))
  ALLOCATE(array2(2))
  array1 = (/"O","OH","HO2"/)
  array2 = (/"O2    ","HO2H    "/)

  IF ( ANY( array1 .EQ. array2 ) ) THEN
    PRINT *, "DING!!"
  ELSE
    PRINT *, "No DING!!"
  ENDIF
END PROGRAM test
