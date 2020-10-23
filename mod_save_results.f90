MODULE save_results
USE global_functions
USE global_variables
IMPLICIT NONE

CONTAINS

SUBROUTINE save_results_shingledecker
IMPLICIT NONE
INTEGER :: i,j
CHARACTER*20 :: abfile, csvfile
CHARACTER*30 :: abform, csvform
CHARACTER*10 :: species
CHARACTER*40 :: keycomment
DOUBLE PRECISION :: abundance

PRINT*, "SAVE_RESULTS: Saving results..."

CALL SYSTEM("rm -rf ab")
CALL SYSTEM("rm -rf csv")
CALL SYSTEM("mkdir ab")
CALL SYSTEM("mkdir csv")
CALL SYSTEM("rm final_gas_abundances.out")
CALL SYSTEM("rm final_surface_abundances.out")
CALL SYSTEM("rm final_bulk_abundances.out")

1004 FORMAT(a10, 14x, e14.12, 1x, a40)
1957 FORMAT(2E10.3)
1958 FORMAT(E10.4,(A1),E10.4)

keycomment = ";"

DO i=1,nspecies
  abfile  = TRIM(s(i)%name)//".ab"
  csvfile = TRIM(s(i)%name)//".csv"
  OPEN ( UNIT=20, FILE=abfile,STATUS='unknown',ACCESS='append' )
  OPEN ( UNIT=30, FILE=csvfile,STATUS='unknown',ACCESS='append' )
  OPEN ( UNIT=40, FILE="final_gas_abundances.out",STATUS='unknown',ACCESS='append' )
  OPEN ( UNIT=50, FILE="final_surface_abundances.out",STATUS='unknown',ACCESS='append' )
  OPEN ( UNIT=60, FILE="final_bulk_abundances.out",STATUS='unknown',ACCESS='append' )
  WRITE(20,'(A103)') "! time [year]; Each column is the abundance (relative to H) [number ratio] for several spatial positions"
  WRITE(30,'(A103)') "! time [year]; Each column is the abundance (relative to H) [number ratio] for several spatial positions"
  IF ( MODEL_EXPERIMENT .EQ. 0 ) THEN
    DO j=1,timesteps
      abundance = ABS(s(i)%abundance_out(j)/gdens*ddens)
      IF (abundance .LT. 1.0e-35) abundance = 0.0E0
      WRITE(20,1957) timesteps_out(j)/year,abundance
      WRITE(30,1958) timesteps_out(j)/year,",",abundance
    ENDDO

    IF (s(i)%name(1:1) .EQ. 'g') THEN
      WRITE(50,1004) s(i)%name,s(i)%abundance_out(j)/gdens*ddens, keycomment
    ELSE IF (s(i)%name(1:1) .EQ. 'b') THEN
      WRITE(60,1004) s(i)%name,s(i)%abundance_out(j)/gdens*ddens, keycomment
    ELSE
      WRITE(40,1004) s(i)%name,s(i)%abundance_out(j)/gdens*ddens, keycomment
    ENDIF
  ELSE
    DO j=1,timesteps
!      abundance = s(i)%abundance_out(j)/s(species_idx('bH2O      '))%abundance_out(j)
!      abundance = s(i)%abundance_out(j)/(2.0*1.0e20*1.0e-12)
      abundance = s(i)%abundance_out(j)
      IF (abundance .LT. 1.0e-35) abundance = 0.0E0
      WRITE(20,1957) timesteps_out(j),abundance
      WRITE(30,1958) timesteps_out(j),",",abundance
    ENDDO
  ENDIF

  CLOSE (20)
  CLOSE (30)
  CLOSE (40)
  CLOSE (50)
  CLOSE (60)
ENDDO

CALL SYSTEM("mv *.ab ab/")
CALL SYSTEM("mv *.csv csv/")

PRINT*, "SAVE_RESULTS: Done"

END SUBROUTINE save_results_shingledecker

SUBROUTINE save_results_semenov
IMPLICIT NONE
integer i,j,l,nait,nlen,li,lt,iana,is,iff,lx
character*5 ait

PRINT*, "SAVE_RESULTS: Saving results..."


!==============================================================================
! WRITE results to the output files:
!==============================================================================

outfile_nameprefix = 'results'
nlen = LEN_TRIM(outfile_nameprefix)
! Open output file '*.idl':

IF (is_disk_model==1) THEN

  IF (curpoint.LE.9) THEN
    WRITE (ait,'(I1)') curpoint
        nait = 1
    ELSE IF (curpoint.LE.99) THEN
        WRITE (ait,'(I2)') curpoint
        nait = 2
    ELSE  IF (curpoint.LE.999) THEN
        WRITE (ait,'(I3)') curpoint
        nait = 3
    ELSE  IF (curpoint.LE.9999) THEN
        WRITE (ait,'(I4)') curpoint
        nait = 4
    ELSE
        WRITE (ait,'(I5)') curpoint
        nait = 5
    END IF

!! Open output file 'out':
!        CALL len_tri2(outfile_nameprefix,30,nlen)

    OPEN (unit=20,file=outfile_nameprefix(1:nlen)//'_'//ait(1:nait)//'.out',status='unknown',access='append')
    REWIND 20
    OPEN (unit=30,file=outfile_nameprefix(1:nlen)//'_'//ait(1:nait)//'.idl',status='unknown',access='append')
    REWIND 30

ELSE

  OPEN (20, FILE=outfile_nameprefix(1:nlen)//'.out', STATUS='UNKNOWN', ACCESS='APPEND')
  REWIND (20)
  OPEN (30, FILE=outfile_nameprefix(1:nlen)//'.idl', STATUS='UNKNOWN', ACCESS='APPEND')
  REWIND (30)

ENDIF

! 1) Write current parameters to the 'out//.idl':
WRITE(30,'(A18)') '# ART input file: '
WRITE(30,'(a80)') chem_file
WRITE(30,'(a80)') !'dummy line'
WRITE(30,'(I5)') curpoint
WRITE(30,'(1PE9.2)') rs
WRITE(30,'(1PE9.2)') zs
WRITE(30,'(0PF5.0)') T
WRITE(30,'(1PE9.2)') rho
WRITE(30,'(1PE9.2)') drho
WRITE(30,'(1PE9.2)') dust2gas
WRITE(30,'(1PE9.2)') agr
WRITE(30,'(1PE9.2)') AvSt
WRITE(30,'(1PE9.2)') AvIS
WRITE(30,'(1PE9.2)') G0_stellar
WRITE(30,'(1PE9.2)') ZetaCR
WRITE(30,'(1PE9.2)') ZetaX
WRITE(30,'(1PE9.2)') albedo_UV
WRITE(30,'(1PE9.2)') tend
WRITE(30,'(1PE9.2)') tstart
WRITE(30,'(I4)') init_non_zero
WRITE(30,'(A10,1x,D22.15)') (s_init(j)%name, s_init(j)%abundance, j = 1, init_non_zero)
WRITE(30,*) nspecies
WRITE(30,'((9A10),:)') (s(j)%name, j = 1, nspecies)
WRITE(30,*) timesteps
WRITE(30,*) (timesteps_out(j), j = 1, timesteps)
WRITE(30,*) nreactions
WRITE(30,'(10D12.5)') (r(j)%rate, j = 1, nreactions)
WRITE(30,'(10D12.5)') ((dmax1(s(j)%abundance_out(l)/gdens*ddens, 1.0D-99), l=1,timesteps), j=1,nspecies)
CLOSE(30)

! 2) Write current parameters to the 'out':
li = nspecies
lt = li+1
iana =1
WRITE(20,20)
WRITE(20,75)
75 FORMAT(21x,39('*'))
WRITE(20,84) lt-1
84 format(20X,' **  CIRCUMSTELLAR ENVELOPE CHEMISTRY **',/, &
   20X,' **          ',1I3,' SPECIES SET          **',/, &
   20X,' **            F77 VERSION            **',/, &
   20X,' **             31/03/2007            **')
write(20,75)
write(20,20)
21 format(1x,' SPECIES CONTAINED IN SCHEME: ',//)
write(20,22) (s(j)%name,j=1,nspecies)
22 format(8(1x,a10))
write(20,23)
23 format(/)
write(20,36) nspecies
36 format(1x,' NUMBER OF VALID SPECIES USED = ',1i4)
write(20,45) nreactions
45 format(1x,' NUMBER OF VALID REACTIONS USED = ',1i6)
write(20,23)
20 format(//)

! Write input parameters:
write(20,64)
64 format(3x,' INITIAL VALUES: '/)
WRITE(20,65) curpoint, rs, zs, gdens, T, rho, Tdust, drho, dust2gas, agr, AvSt, AvIS, G0_stellar, (ZetaCR+ZetaX), albedo_UV, tstart, tend, 0.0d0
65 FORMAT( &
   3X,' Grid point  = ',I5,/, &
   3X,' Radius      = ',1PE11.3,' AU',/, &
   3X,' Height      = ',1PE11.3,' AU',/, &
   3X,' n(H+2H2)    = ',1PE11.3,' cm^(-3)',/, &
   3X,' Tg          = ',0PF8.1,'K',/, &
   3X,' rho_g       = ',1PE11.3,' g/cm^3',/, &
   3X,' Td          = ',0PF8.1,'K',/, &
   3X,' rho_d       = ',1PE11.3,' g/cm^3',/, &
   3X,' Mdust/Mgas  = ',1PE11.3,/, &
   3X,' grain size  = ',1PE11.3,' cm',/, &
   3X,' Av(stellar) = ',1PE11.3,' mag.',/, &
   3X,' Av(IS)      = ',1PE11.3,' mag.',/, &
   3X,' G0(stellar) = ',1PE11.3,' G0(IS)',/, &
   3X,' Zeta        = ',1PE11.3,/, &
   3X,' albedo(UV)  = ',1PE11.3,/, &
   3X,' Start time  = ',0PF9.0,' years',/, &
   3X,' Finish time = ',0PF9.0,' years',/, &
   3X,' Diff. coeff.= ',1PE11.3)
WRITE(20,20)
! Write calculated abundances:
write(20,63)
63 format(3X,' CALCULATED ABUNDANCES: '/)
is = 1
iff = min(6,nspecies)
32 write(20,41) (s(j)%name,j=is,iff)
write(20,76)
do l = 1, timesteps
  write(20,30) timesteps_out(l)/year,(dmax1(s(j)%abundance_out(l)/gdens*ddens, 1.0D-99),j=is,iff)
end do
30 FORMAT(1x,7(1pe11.3))
WRITE(20,41)(s(j)%name,j=is,iff)
41 FORMAT(6x,'time',6x,6(1a10,1x))
write(20,20)
is=is+6
iff=iff+6
! Interupt?
if (iff.gt.nspecies) iff=nspecies
if (is.ge.li) go to 38
go to 32
! Last output statement:
38 lx=li
lt=lt-1
76 FORMAT(1x,80('-'))
CLOSE (20)

END SUBROUTINE save_results_semenov

SUBROUTINE save_results_bulk
IMPLICIT NONE
INTEGER :: i, j, l ,nait,nlen,li,lt,iana,is,iff,lx, nml_max_n, nml_n
CHARACTER*3 :: nml_max_s, nml_s

nml_max_n = aint(nml_max)
IF (nml_max<1000) WRITE(nml_max_s,'(i3)')nml_max_n
IF (nml_max<100) WRITE(nml_max_s,'(i2)')nml_max_n
IF (nml_max<10) WRITE(nml_max_s,'(i1)')nml_max_n

nml_n = aint(nml)
IF (nml<1000) WRITE(nml_s,'(i3)')nml_n
IF (nml<100) WRITE(nml_s,'(i2)')nml_n
IF (nml<10) WRITE(nml_s,'(i1)')nml_n

OPEN (200, FILE='results_bulk.'//nml_max_s(1:LEN_TRIM(nml_max_s))//'.'//nml_s(1:LEN_TRIM(nml_s))//'.out', STATUS='UNKNOWN', ACCESS='APPEND')
REWIND (200)
OPEN (300, FILE='results_bulk.'//nml_max_s(1:LEN_TRIM(nml_max_s))//'.'//nml_s(1:LEN_TRIM(nml_s))//'.idl', STATUS='UNKNOWN', ACCESS='APPEND')
REWIND (300)

! 1) Write current parameters to the 'out//.idl':
WRITE(300,'(A18)') '# ART input file: '
WRITE(300,'(a80)') chem_file
WRITE(300,'(a80)') !'dummy line'
WRITE(300,'(I5)') curpoint
WRITE(300,'(1PE9.2)') rs
WRITE(300,'(1PE9.2)') zs
WRITE(300,'(0PF5.0)') T
WRITE(300,'(1PE9.2)') rho
WRITE(300,'(1PE9.2)') drho
WRITE(300,'(1PE9.2)') dust2gas
WRITE(300,'(1PE9.2)') agr
WRITE(300,'(1PE9.2)') AvSt
WRITE(300,'(1PE9.2)') AvIS
WRITE(300,'(1PE9.2)') G0_stellar
WRITE(300,'(1PE9.2)') ZetaCR
WRITE(300,'(1PE9.2)') ZetaX
WRITE(300,'(1PE9.2)') albedo_UV
WRITE(300,'(1PE9.2)') tend
WRITE(300,'(1PE9.2)') tstart
WRITE(300,'(I4)') init_non_zero
WRITE(300,'(A10,1x,D22.15)') (s_init(j)%name, s_init(j)%abundance, j = 1, init_non_zero)
WRITE(300,*) nspecies
WRITE(300,'((9A10),:)') (s(j)%name, j = 1, nspecies)
WRITE(300,*) nml, nml_max
!WRITE(300,*) (timesteps_nml(j), j = 1, nml_max)
WRITE(300,*) nreactions
WRITE(300,'(10D12.5)') (r(j)%rate, j = 1, nreactions)
WRITE(300,'(10D17.10)') ((abundances_bulk(j,l)/gdens*ddens, l=1,nml_max), j=1,nspecies)
CLOSE (300)

! 2) Write current parameters to the 'out':
li = nspecies
lt = li+1
iana =1
WRITE(200,20)
WRITE(200,75)
75 FORMAT(21x,39('*'))
WRITE(200,84) lt-1
84 format(20X,' **  CIRCUMSTELLAR ENVELOPE CHEMISTRY **',/, &
   20X,' **          ',1I3,' SPECIES SET          **',/, &
   20X,' **            F77 VERSION            **',/, &
   20X,' **             31/03/2007            **')
write(200,75)
write(200,20)
21 format(1x,' SPECIES CONTAINED IN SCHEME: ',//)
write(200,22) (s(j)%name,j=1,nspecies)
22 format(8(1x,a10))
write(200,23)
23 format(/)
write(200,36) nspecies
36 format(1x,' NUMBER OF VALID SPECIES USED = ',1i4)
write(200,45) nreactions
45 format(1x,' NUMBER OF VALID REACTIONS USED = ',1i6)
write(200,23)
20 format(//)

! Write input parameters:
write(200,64)
64 format(3x,' INITIAL VALUES: '/)
WRITE(200,65) curpoint, rs, zs, gdens, T, rho, Tdust, drho, dust2gas, agr, AvSt, AvIS, G0_stellar, (ZetaCR+ZetaX), albedo_UV, tstart, tend, 0.0d0
65 FORMAT( &
   3X,' Grid point  = ',I5,/, &
   3X,' Radius      = ',1PE11.3,' AU',/, &
   3X,' Height      = ',1PE11.3,' AU',/, &
   3X,' n(H+2H2)    = ',1PE11.3,' cm^(-3)',/, &
   3X,' Tg          = ',0PF8.1,'K',/, &
   3X,' rho_g       = ',1PE11.3,' g/cm^3',/, &
   3X,' Td          = ',0PF8.1,'K',/, &
   3X,' rho_d       = ',1PE11.3,' g/cm^3',/, &
   3X,' Mdust/Mgas  = ',1PE11.3,/, &
   3X,' grain size  = ',1PE11.3,' cm',/, &
   3X,' Av(stellar) = ',1PE11.3,' mag.',/, &
   3X,' Av(IS)      = ',1PE11.3,' mag.',/, &
   3X,' G0(stellar) = ',1PE11.3,' G0(IS)',/, &
   3X,' Zeta        = ',1PE11.3,/, &
   3X,' albedo(UV)  = ',1PE11.3,/, &
   3X,' Start time  = ',0PF9.0,' years',/, &
   3X,' Finish time = ',0PF9.0,' years',/, &
   3X,' Diff. coeff.= ',1PE11.3)
WRITE(200,20)
! Write calculated abundances:
write(200,63)
63 format(3X,' CALCULATED ABUNDANCES: '/)
is = 1
iff = min(6,nspecies)
32 write(200,41) (s(j)%name,j=is,iff)
write(200,76)
do l = 1, nml_max
!  write(200,30) timesteps_nml(l)/year,(abundances_bulk(j,l)/gdens*ddens,j=is,iff)
end do
30 FORMAT(1x,7(1pe11.3))
WRITE(200,41)(s(j)%name,j=is,iff)
41 FORMAT(6x,'time',6x,6(1a10,1x))
write(200,20)
is=is+6
iff=iff+6
! Interupt?
if (iff.gt.nspecies) iff=nspecies
if (is.ge.li) go to 38
go to 32
! Last output statement:
38 lx=li
lt=lt-1
76 FORMAT(1x,80('-'))

CLOSE (200)

END SUBROUTINE save_results_bulk

END
