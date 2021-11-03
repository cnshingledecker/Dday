MODULE read_rate06
USE global_variables
USE global_functions
IMPLICIT NONE

CONTAINS

!Subroutine to read chemical database file
SUBROUTINE read_rate06database
IMPLICIT NONE
INTEGER :: i, ii, j, jj, idx, ir1, ir2, ip1, ip2, ip3, ip4, ip5, dummy, s_r_counter
INTEGER :: first_suprathermal_react,first_suprathermal_species
INTEGER :: Nsup_g,Nlines,NGRIONS
INTEGER :: io
INTEGER :: prodatoms,reactatoms
INTEGER :: NTOTALSPECIES,NSUPRATHERMAL
INTEGER :: first_surface_suprathermal
REAL*8  :: a, b, c
REAL*8  :: apriori_nml
CHARACTER*10 :: groundstate,s_name, r1, r2, p1, p2, p3, p4, p5
TYPE (reaction), DIMENSION(:), ALLOCATABLE :: rtemp

! Set a priori number of monolayers
apriori_nml = 1000000 !10.0*(ICE_THICK/5.0e-8)

!Pre-reading of ratefile to find the numbers of surface species and surface reactions
n_surf_spec = 0
n_surf_react = 0

OPEN(1, FILE=chem_file, STATUS='OLD', ERR=100)
READ(1,*)nspecies

NGRIONS = 0
DO i = 1, nspecies
  READ(1,'(a10)')s_name
  IF (s_name(1:1) == 'g') n_surf_spec = n_surf_spec + 1

  IF ( s_name(1:1)=='g' ) THEN
    IF ( s_name(LEN_TRIM(s_name):LEN_TRIM(s_name)) .EQ. '+' .OR. s_name(LEN_TRIM(s_name):LEN_TRIM(s_name)) .EQ. '-') THEN
      NGRIONS = NGRIONS + 1
    ENDIF
  ENDIF


  PRINT *, "The species is ",s_name
ENDDO

! Calculate number of suprathermal surface species
NSUPRATHERMAL = n_surf_spec - NGRIONS
PRINT *, "The number of suprathermal species is ", NSUPRATHERMAL

!Find the number of charged grain species
DO i = 1, nspecies
ENDDO

PRINT *, 'NGRIONS =',NGRIONS


READ(1,*)nreactions

!Counting reactions with 'g'-species, to be able to add corresponding bulk reactions
DO i = 1, nreactions
    READ(1,1000)idx, r1, r2, p1, p2, p3, p4, p5, a, b, c
    IF (r1(1:1)=='g' .AND. r2/='FREEZE' .AND. r2/='DESORB' .AND. p1(1:1)=='g') n_surf_react = n_surf_react + 1
ENDDO

CLOSE (1)

!Main reading of ratefile
OPEN(1, FILE=chem_file, STATUS='OLD', ERR=100)


! C. N. Shingledecker
IF ( SUPRATHERMAL .EQ. 1) THEN
! 1x for the suprathermal surface species - Number of charged species
! 1x for the suprathermal bulk species - Number of charged species
! In total: n_grain_species + 2*(n_grain_species - NGRIONS)
  NTOTALSPECIES = nspecies + 3*n_surf_spec - 2*NGRIONS
ELSE
! 1x for the bulk species
  NTOTALSPECIES = nspecies + n_surf_spec
ENDIF



READ(1,*)dummy
IF ( suprathermal .EQ. 1 ) THEN
  ALLOCATE ( s(NTOTALSPECIES) )
  IF ( MODEL_EXPERIMENT .EQ. 0 ) THEN
    ALLOCATE (abundances_bulk(NTOTALSPECIES,10000)) !Up to 10000 monolayers
  ELSE
    ALLOCATE (abundances_bulk(NTOTALSPECIES,INT(apriori_nml))) !Arbitrarily many monolayers
  ENDIF
  DO i = 1,NTOTALSPECIES 
     ALLOCATE (s(i)%abundance_out(timesteps))
     s(i)%abundance_out(:) = 0.0d0
  ENDDO
ELSE
  ALLOCATE (s(NTOTALSPECIES))
  IF ( MODEL_EXPERIMENT .EQ. 0 ) THEN
    ALLOCATE (abundances_bulk(NTOTALSPECIES,500)) !Up to 300 monolayers
  ELSE
    ALLOCATE (abundances_bulk(NTOTALSPECIES,INT(apriori_nml))) !Arbitrarily many monolayers
  ENDIF
  DO i = 1,NTOTALSPECIES
     ALLOCATE (s(i)%abundance_out(timesteps))
     s(i)%abundance_out(:) = 0.0d0
  ENDDO
END IF
IF ( MODEL_EXPERIMENT .EQ. 0 ) THEN
  ALLOCATE (timesteps_nml(3000))
ELSE
  ALLOCATE (timesteps_nml(INT(apriori_nml)))
ENDIF
s(:)%edes           = 0.0d0
s(:)%racc           = 0.0d0
s(:)%rdes           = 0.0d0
s(:)%abundance      = 0.0d0
s(:)%frac_abundance = 0.0d0


first_surf_spec = 0

DO i = 1, nspecies
  READ(1,'(a10)')s(i)%name
  PRINT *, s(i)%name
  s(i)%idx = i
  s(i)%gas_idx = i
  s(i)%weight = aweight(s(i)%name)
  PRINT *, "Weight = ", s(i)%weight
  s(i)%natoms = numatoms(s(i)%name)
  PRINT *, "N_atoms = ", s(i)%natoms
  IF (first_surf_spec == 0 .AND. s(i)%name(1:1) == 'g') first_surf_spec = i
ENDDO



!Looking for indexes of gas counterparts of surface species
DO i = 1, nspecies
  IF (s(i)%name(1:1)=='g') THEN
    DO j = 1, nspecies
      IF (s(j)%name==s(i)%name(2:LEN_TRIM(s(i)%name))) s(i)%gas_idx = j
      s(i)%enthalpia = 0.0d0
      s(i)%enthalpia_known = 0
    ENDDO
  ENDIF
ENDDO

DO i = 1, nspecies
  IF (s(i)%name(1:1)=='g') THEN
    PRINT *, s(i)%name, " has gas counterpart ", s(s(i)%gas_idx)%name
  ENDIF
ENDDO




IF ( suprathermal .EQ. 1 ) THEN
  j = 1
  PRINT *, "Now adding surface suprathermal species"
  PRINT *, "first_surf_spec =", s(first_surf_spec)%name, " at idx ", first_surf_spec
  PRINT *, "There are currently nspecies= ",nspecies," with the last being ", s(nspecies)%name
  PRINT *, "The size of the s-array is ",SIZE(s)
  !Adding suprathermal surface species
  DO i = first_surf_spec, nspecies
    PRINT *, "In loop"
    PRINT *, s(i)%name, " at idx = ",i
    IF ( s(i)%name(LEN_TRIM(s(i)%name):LEN_TRIM(s(i)%name)) .NE. '+' .AND. &
         s(i)%name(LEN_TRIM(s(i)%name):LEN_TRIM(s(i)%name)) .NE. '-') THEN
      PRINT *, "Should be made suprathermal"
      PRINT *, "j =",j," nspecies+j=",nspecies+j," and should be GE ",nspecies+1
      s(nspecies+j)%name = s(i)%name(1:LEN_TRIM(s(i)%name))//'*'
      s(nspecies+j)%idx = nspecies + j
      s(nspecies+j)%gas_idx = s(i)%gas_idx
      s(nspecies+j)%weight = s(i)%weight
      s(nspecies+j)%natoms = s(i)%natoms
      s(nspecies+j)%enthalpia_known = s(i)%enthalpia_known
      s(nspecies+j)%enthalpia = s(i)%enthalpia
      s(nspecies+j)%edes = s(i)%edes
      s(nspecies+j)%racc = s(i)%racc
      PRINT *, TRIM(s(nspecies+j)%name), " with idx ", s(nspecies+j)%idx, " is the suprathermal counterpart of ",s(i)%name, " at idx ", i
      j = j + 1
    ELSE
      PRINT *, "Should not be made suprathermal"
    ENDIF 
    PRINT *, "*************************************"
  ENDDO
  first_surface_suprathermal = nspecies + 1
  nspecies = nspecies + n_surf_spec - NGRIONS
END IF


PRINT *, "End of adding surface suprathermal species"
PRINT *, s(first_surface_suprathermal)%name


CALL read_enthalpias

READ(1,*)nreactions
IF (bulk_chemistry == 0) THEN
    ALLOCATE (r(nreactions))
ELSE
    ALLOCATE (r(nreactions+n_surf_react))
ENDIF
ALLOCATE (mre_terms(nreactions+n_surf_react))
ALLOCATE (rd_v2_terms(5,nreactions+n_surf_react))
rd_v2_terms(:,:) = 0.0d0
ii = 1
!new_nreactions = nreactions
first_surfreact = 0
Nsup_g          = 0

DO i = 1, nreactions
  READ(1,1000)r(ii)%idx, r(ii)%r1, r(ii)%r2, r(ii)%p1, r(ii)%p2, r(ii)%p3, r(ii)%p4, r(ii)%p5, r(ii)%alpha, r(ii)%beta, r(ii)%gamma
  r(ii)%ir1 = species_idx(r(ii)%r1)
  r(ii)%ir2 = species_idx(r(ii)%r2)
  r(ii)%ip1 = species_idx(r(ii)%p1)
  r(ii)%ip2 = species_idx(r(ii)%p2)
  r(ii)%ip3 = species_idx(r(ii)%p3)
  r(ii)%ip4 = species_idx(r(ii)%p4)
  r(ii)%ip5 = species_idx(r(ii)%p5)
  IF (r(ii)%ir1==-1 .OR. r(ii)%ir2==-1 .OR. r(ii)%ip1==-1 .OR. r(ii)%ip2==-1 .OR. r(ii)%ip3==-1 .OR. r(ii)%ip4==-1 .OR. r(ii)%ip5==-1) THEN
      PRINT '(a31, i5, 1x, 7a10)','Reaction with unknown species: ', r(ii)%idx, r(ii)%r1, r(ii)%r2, r(ii)%p1, r(ii)%p2, r(ii)%p3, r(ii)%p4, r(ii)%p5
      PRINT *, "Around line 150"

      IF (r(ii)%ir1==-1) THEN
	PRINT *, "Problem is with ", r(ii)%r1
      ENDIF

      IF (r(ii)%ir2==-1) THEN
	PRINT *, "Problem is with ", r(ii)%r2
      ENDIF

      IF (r(ii)%ip1==-1) THEN
	PRINT *, "Problem is with ", r(ii)%p1
	PRINT *, "Weight =", s(r(ii)%ip1)%weight
	PRINT *, "Num. atoms =", s(r(ii)%ip1)%natoms
	PRINT *, "Idx =", s(r(ii)%ip1)%idx
	PRINT *, "Gas Idx =", s(r(ii)%ip1)%gas_idx
      ENDIF

      IF (r(ii)%ip2==-1) THEN
	PRINT *, "Problem is with ", r(ii)%p2
      ENDIF

      IF (r(ii)%ip3==-1) THEN
	PRINT *, "Problem is with ", r(ii)%p3
      ENDIF

      IF (r(ii)%ip4==-1) THEN
	PRINT *, "Problem is with ", r(ii)%p4
      ENDIF

      IF (r(ii)%ip5==-1) THEN
	PRINT *, "Problem is with ", r(ii)%p5
      ENDIF

      STOP
  ENDIF

  CALL get_reaction_thermodynamics(ii)

  r(ii)%rtype = get_rtype(r(ii)%r1,r(ii)%r2)
  IF (r(ii)%rtype == 12) s(r(ii)%ir1)%edes = r(ii)%gamma
  IF (r(ii)%rtype == 12) s(s(r(ii)%ir1)%gas_idx)%edes = r(ii)%gamma

  IF (first_surfreact == 0 .AND. r(ii)%rtype == 13) first_surfreact = ii

    CALL get_rd_efficiency(ii)


    write(37,1001)r(ii)%idx,r(ii)%rtype,r(ii)%r1,r(ii)%r2,r(ii)%p1,r(ii)%p2,r(ii)%p3,r(ii)%p4,r(ii)%p5,r(ii)%exothermicity,r(ii)%exothermicity_known, r(ii)%alpha !(des_reactive*P)/(1+des_reactive*P)
  !Assigning the desorption energies to the surface species


  ii = ii + 1
ENDDO

CLOSE (1)



!Adding bulk species
first_bulk_spec = nspecies + 1
DO i = 1, n_surf_spec
    s(nspecies+i)%name = 'b'//s(first_surf_spec+i-1)%name(2:LEN_TRIM(s(first_surf_spec+i-1)%name))
    s(nspecies+i)%idx = nspecies + i
    s(nspecies+i)%gas_idx = s(first_surf_spec+i-1)%gas_idx
    s(nspecies+i)%weight = s(first_surf_spec+i-1)%weight
    s(nspecies+i)%natoms = s(first_surf_spec+i-1)%natoms
    s(nspecies+i)%enthalpia_known = s(first_surf_spec+i-1)%enthalpia_known
    s(nspecies+i)%enthalpia = s(first_surf_spec+i-1)%enthalpia
    s(nspecies+i)%edes = s(first_surf_spec+i-1)%edes
    s(nspecies+i)%racc = s(first_surf_spec+i-1)%racc
ENDDO
nspecies = nspecies + n_surf_spec

first_suprathermal_species = nspecies+1

PRINT *, "First bulk species at ",first_bulk_spec," = ",TRIM(s(first_bulk_spec)%name)

! C. N. Shingledecker
IF ( suprathermal .EQ. 1 ) THEN
  !Adding suprathermal bulk species
  j = 1
  DO i = 1, n_surf_spec
    IF ( s(first_surf_spec+i-1)%name(LEN_TRIM(s(first_surf_spec+i-1)%name):LEN_TRIM(s(first_surf_spec+i-1)%name)) .NE. '+' .AND. &
         s(first_surf_spec+i-1)%name(LEN_TRIM(s(first_surf_spec+i-1)%name):LEN_TRIM(s(first_surf_spec+i-1)%name)) .NE. '-') THEN
      s(nspecies+j)%name = 'b'//s(first_surf_spec+i-1)%name(2:LEN_TRIM(s(first_surf_spec+i-1)%name))//'*'
      s(nspecies+j)%idx = nspecies + i
      s(nspecies+j)%gas_idx = s(first_surf_spec+i-1)%gas_idx
      s(nspecies+j)%weight = s(first_surf_spec+i-1)%weight
      s(nspecies+j)%natoms = s(first_surf_spec+i-1)%natoms
      s(nspecies+j)%enthalpia_known = s(first_surf_spec+i-1)%enthalpia_known
      s(nspecies+j)%enthalpia = s(first_surf_spec+i-1)%enthalpia
      s(nspecies+j)%edes = s(first_surf_spec+i-1)%edes
      s(nspecies+j)%racc = s(first_surf_spec+i-1)%racc
      PRINT *, s(nspecies+j)%name, "j = ",j, "nspecies+j=",nspecies+j
      j = j + 1
    ELSE
      PRINT *, s(first_surf_spec+i-1)%name, " at ", first_surf_spec+i-1, " should not be made suprathermal"
    ENDIF 
  ENDDO
  nspecies = nspecies + n_surf_spec - NGRIONS

  DO i = 1, nspecies
    PRINT *, s(i)%name, " at i = ", i
  END DO



END IF

DO i = 1, nspecies
  IF (s(i)%name(1:1)=='g') THEN
    PRINT *, s(i)%name, " has gas counterpart ", s(s(i)%gas_idx)%name
    PRINT *, "And bulk counterpart ",s(get_bulk_idx(s(i)%name))%name
    PRINT *, "******************************************"
  ENDIF
ENDDO



!Adding bulk reactions
first_bulkreact = nreactions + 1
IF (bulk_chemistry>0) THEN
    s_r_counter = 0
    DO i = 1, nreactions
        IF (r(i)%r1(1:1)=='g' .AND. r(i)%r2/='FREEZE' .AND. r(i)%r2/='DESORB' .AND. r(i)%p1(1:1)=='g') THEN
            s_r_counter = s_r_counter + 1
            r(nreactions+s_r_counter)%idx = nreactions+s_r_counter

            r(nreactions+s_r_counter)%r1 = r(i)%r1
            r(nreactions+s_r_counter)%r2 = r(i)%r2
            r(nreactions+s_r_counter)%p1 = r(i)%p1
            r(nreactions+s_r_counter)%p2 = r(i)%p2
            r(nreactions+s_r_counter)%p3 = r(i)%p3
            r(nreactions+s_r_counter)%p4 = r(i)%p4
            r(nreactions+s_r_counter)%p5 = r(i)%p5

            IF (r(i)%r1(1:1)=='g') r(nreactions+s_r_counter)%r1 = 'b'//r(i)%r1(2:LEN_TRIM(r(i)%r1))
            IF (r(i)%r2(1:1)=='g') r(nreactions+s_r_counter)%r2 = 'b'//r(i)%r2(2:LEN_TRIM(r(i)%r2))
            IF (r(i)%p1(1:1)=='g') r(nreactions+s_r_counter)%p1 = 'b'//r(i)%p1(2:LEN_TRIM(r(i)%p1))
            IF (r(i)%p2(1:1)=='g') r(nreactions+s_r_counter)%p2 = 'b'//r(i)%p2(2:LEN_TRIM(r(i)%p2))
            IF (r(i)%p3(1:1)=='g') r(nreactions+s_r_counter)%p3 = 'b'//r(i)%p3(2:LEN_TRIM(r(i)%p3))
            IF (r(i)%p4(1:1)=='g') r(nreactions+s_r_counter)%p4 = 'b'//r(i)%p4(2:LEN_TRIM(r(i)%p4))
            IF (r(i)%p5(1:1)=='g') r(nreactions+s_r_counter)%p5 = 'b'//r(i)%p5(2:LEN_TRIM(r(i)%p5))

            r(nreactions+s_r_counter)%ir1 = species_idx(r(nreactions+s_r_counter)%r1)
            r(nreactions+s_r_counter)%ir2 = species_idx(r(nreactions+s_r_counter)%r2)
            r(nreactions+s_r_counter)%ip1 = species_idx(r(nreactions+s_r_counter)%p1)
            r(nreactions+s_r_counter)%ip2 = species_idx(r(nreactions+s_r_counter)%p2)
            r(nreactions+s_r_counter)%ip3 = species_idx(r(nreactions+s_r_counter)%p3)
            r(nreactions+s_r_counter)%ip4 = species_idx(r(nreactions+s_r_counter)%p4)
            r(nreactions+s_r_counter)%ip5 = species_idx(r(nreactions+s_r_counter)%p5)

            r(nreactions+s_r_counter)%alpha = r(i)%alpha
            r(nreactions+s_r_counter)%beta  = r(i)%beta
            r(nreactions+s_r_counter)%gamma = r(i)%gamma

            r(nreactions+s_r_counter)%rtype = get_rtype(r(nreactions+s_r_counter)%r1,r(nreactions+s_r_counter)%r2)

            r(nreactions+s_r_counter)%exothermicity_known = r(i)%exothermicity_known
            r(nreactions+s_r_counter)%exothermicity = r(i)%exothermicity
        ENDIF
    ENDDO

    nreactions = nreactions + n_surf_react

    DO i = 1, nreactions
        write(38,1001)r(i)%idx,r(i)%rtype,r(i)%r1,r(i)%r2,r(i)%p1,r(i)%p2,r(i)%p3,r(i)%p4,r(i)%p5,r(i)%exothermicity,r(i)%exothermicity_known, r(i)%alpha !(des_reactive*P)/(1+des_reactive*P)
    ENDDO

ENDIF

! Count how many suprathermal reactions we have to add
Nsup_g = 0
DO i = 1,nreactions
  IF (r(i)%rtype .EQ. 13 .OR. r(i)%rtype .EQ. 14) THEN
    ! If both reactants are the same, just add one suprathermal
    ! reaction. Otherwise, add two, one for each reactant.
    IF ( r(i)%r1 .EQ. r(i)%r2 ) THEN
      Nsup_g = Nsup_g + 1
    ELSE
      Nsup_g = Nsup_g + 2
    END IF
  ENDIF
ENDDO

! Now generate new reactions array with additional space for suprathermal
! reactions.
PRINT *, "original size(r)= ",size(r)
ALLOCATE( rtemp(SIZE(r) + Nsup_g) )
rtemp(1:SIZE(r)) = r
DEALLOCATE( r )
ALLOCATE( r(SIZE(rtemp)) )
r = rtemp
DEALLOCATE( rtemp )

PRINT *, "Nsup_g = ",Nsup_g
PRINT *, "next size(r)= ",size(r)

! C. N. Shingledecker
! Add suprathermal versions of thermal surface and bulk reactions involving 1 suprathermal reactant
IF ( suprathermal .EQ. 1 ) THEN
    s_r_counter = 0
    first_suprathermal_react = nreactions + 1
    DO i = first_surfreact, nreactions
        IF (r(i)%rtype .EQ. 13 .OR. r(i)%rtype .EQ. 14) THEN

          ! Just make one copy with suprathermal species
          IF ( r(i)%r1 .EQ. r(i)%r2 ) THEN
            s_r_counter = s_r_counter + 1
            r(nreactions+s_r_counter)%idx = nreactions+s_r_counter

            r(nreactions+s_r_counter)%r1 = r(i)%r1
            r(nreactions+s_r_counter)%r2 = r(i)%r2
            r(nreactions+s_r_counter)%p1 = r(i)%p1
            r(nreactions+s_r_counter)%p2 = r(i)%p2
            r(nreactions+s_r_counter)%p3 = r(i)%p3
            r(nreactions+s_r_counter)%p4 = r(i)%p4
            r(nreactions+s_r_counter)%p5 = r(i)%p5

            ! Make the first reactant suprathermal
            r(nreactions+s_r_counter)%r1 = r(i)%r1(1:LEN_TRIM(r(i)%r1))//'*'

            r(nreactions+s_r_counter)%ir1 = species_idx(r(nreactions+s_r_counter)%r1)
            r(nreactions+s_r_counter)%ir2 = species_idx(r(nreactions+s_r_counter)%r2)
            r(nreactions+s_r_counter)%ip1 = species_idx(r(nreactions+s_r_counter)%p1)
            r(nreactions+s_r_counter)%ip2 = species_idx(r(nreactions+s_r_counter)%p2)
            r(nreactions+s_r_counter)%ip3 = species_idx(r(nreactions+s_r_counter)%p3)
            r(nreactions+s_r_counter)%ip4 = species_idx(r(nreactions+s_r_counter)%p4)
            r(nreactions+s_r_counter)%ip5 = species_idx(r(nreactions+s_r_counter)%p5)

            r(nreactions+s_r_counter)%alpha = r(i)%alpha
            r(nreactions+s_r_counter)%beta  = r(i)%beta
            r(nreactions+s_r_counter)%gamma = r(i)%gamma

            r(nreactions+s_r_counter)%rtype = get_rtype(r(nreactions+s_r_counter)%r1,r(nreactions+s_r_counter)%r2)

            r(nreactions+s_r_counter)%exothermicity_known = 0
            r(nreactions+s_r_counter)%exothermicity = 0.0d0
          ELSE
          ! Make two copies, one with r1* and another with r2*
            DO ii = 1,2
              s_r_counter = s_r_counter + 1
              r(nreactions+s_r_counter)%idx = nreactions+s_r_counter

              PRINT *, "nreactions=",nreactions
              PRINT *, "s_r_counter=",s_r_counter
              PRINT *, "n_surf_react=",n_surf_react
              PRINT *, "nreactions+s_r_counter=",nreactions+s_r_counter
              PRINT *, "size(r)=",size(r,1)
              PRINT *, r(i)%r1," + ",r(i)%r2," -> ",r(i)%p1," + ",r(i)%p2," + ",r(i)%p3," + ",r(i)%p4
              PRINT *, "************************"


              r(nreactions+s_r_counter)%r1 = r(i)%r1
              r(nreactions+s_r_counter)%r2 = r(i)%r2
              r(nreactions+s_r_counter)%p1 = r(i)%p1
              r(nreactions+s_r_counter)%p2 = r(i)%p2
              r(nreactions+s_r_counter)%p3 = r(i)%p3
              r(nreactions+s_r_counter)%p4 = r(i)%p4
              r(nreactions+s_r_counter)%p5 = r(i)%p5

              IF (ii .EQ. 1) r(nreactions+s_r_counter)%r1 = r(i)%r1(1:LEN_TRIM(r(i)%r1))//'*'
              IF (ii .EQ. 2) r(nreactions+s_r_counter)%r2 = r(i)%r2(1:LEN_TRIM(r(i)%r2))//'*'

              r(nreactions+s_r_counter)%ir1 = species_idx(r(nreactions+s_r_counter)%r1)
              r(nreactions+s_r_counter)%ir2 = species_idx(r(nreactions+s_r_counter)%r2)
              r(nreactions+s_r_counter)%ip1 = species_idx(r(nreactions+s_r_counter)%p1)
              r(nreactions+s_r_counter)%ip2 = species_idx(r(nreactions+s_r_counter)%p2)
              r(nreactions+s_r_counter)%ip3 = species_idx(r(nreactions+s_r_counter)%p3)
              r(nreactions+s_r_counter)%ip4 = species_idx(r(nreactions+s_r_counter)%p4)
              r(nreactions+s_r_counter)%ip5 = species_idx(r(nreactions+s_r_counter)%p5)

              r(nreactions+s_r_counter)%alpha = r(i)%alpha
              r(nreactions+s_r_counter)%beta  = r(i)%beta
              r(nreactions+s_r_counter)%gamma = r(i)%gamma

              r(nreactions+s_r_counter)%rtype = get_rtype(r(nreactions+s_r_counter)%r1,r(nreactions+s_r_counter)%r2)

              r(nreactions+s_r_counter)%exothermicity_known = 0
              r(nreactions+s_r_counter)%exothermicity = 0.0d0
            END DO
          END IF
        ENDIF
    ENDDO

    nreactions = nreactions + Nsup_g

    OPEN(1001,FILE="suprathermal_reactions.out",STATUS='REPLACE')
    DO i = first_suprathermal_react, nreactions
        write(1001,1001)r(i)%idx,r(i)%rtype,r(i)%r1,r(i)%r2,r(i)%p1,r(i)%p2,r(i)%p3,r(i)%p4,r(i)%p5,r(i)%exothermicity,r(i)%exothermicity_known, r(i)%alpha !(des_reactive*P)/(1+des_reactive*P)
    ENDDO
    CLOSE(1001)


    ! Now read file with radiolysis reactions
    ! and count the number of processes
    io     = 0
    Nlines = 0
    OPEN(3,FILE='radiolysis.dat',STATUS='OLD',IOSTAT=io)
    DO
      READ(3,*,IOSTAT=io)
      IF (io .NE. 0) EXIT
      Nlines = Nlines + 1
    ENDDO
    CLOSE(3)

    ! Allocate temp reactions object to hold the radiolysis processes
    ALLOCATE( rtemp(nreactions + Nlines) )
    rtemp(1:SIZE(r)) = r

    ! Add the radiolysis processes to the temp reactions object
    OPEN(3,FILE='radiolysis.dat',STATUS='OLD',IOSTAT=io)
    DO ii=nreactions+1,nreactions + Nlines
      READ(3,1000)rtemp(ii)%idx, rtemp(ii)%r1, rtemp(ii)%r2, rtemp(ii)%p1, rtemp(ii)%p2, rtemp(ii)%p3, rtemp(ii)%p4, rtemp(ii)%p5, rtemp(ii)%alpha, rtemp(ii)%beta, rtemp(ii)%gamma
      rtemp(ii)%idx = ii
      rtemp(ii)%ir1 = species_idx(rtemp(ii)%r1)
      rtemp(ii)%ir2 = species_idx(rtemp(ii)%r2)
      rtemp(ii)%ip1 = species_idx(rtemp(ii)%p1)
      rtemp(ii)%ip2 = species_idx(rtemp(ii)%p2)
      rtemp(ii)%ip3 = species_idx(rtemp(ii)%p3)
      rtemp(ii)%ip4 = species_idx(rtemp(ii)%p4)
      rtemp(ii)%ip5 = species_idx(rtemp(ii)%p5)
      IF (rtemp(ii)%ir1==-1 .OR. rtemp(ii)%ir2==-1 .OR. rtemp(ii)%ip1==-1 .OR. rtemp(ii)%ip2==-1 .OR. rtemp(ii)%ip3==-1 .OR. rtemp(ii)%ip4==-1 .OR. rtemp(ii)%ip5==-1) THEN
          PRINT '(a31, i5, 1x, 7a10)','Reaction with unknown species: ', rtemp(ii)%idx, rtemp(ii)%r1, rtemp(ii)%r2, rtemp(ii)%p1, rtemp(ii)%p2, rtemp(ii)%p3, rtemp(ii)%p4, rtemp(ii)%p5
          STOP
      ENDIF

      rtemp(ii)%rtype = get_rtype(rtemp(ii)%r1,rtemp(ii)%r2)
      rtemp(ii)%exothermicity_known = 0.0
      rtemp(ii)%exothermicity = 0.0
    END DO
    CLOSE(3)

    ! Print out radiolysis reactions
    OPEN(1001,FILE="radiolysis_reactions.out",STATUS='REPLACE')
    DO i=nreactions+1,nreactions + Nlines
        write(1001,1001)rtemp(i)%idx,rtemp(i)%rtype,rtemp(i)%r1,rtemp(i)%r2,rtemp(i)%p1,rtemp(i)%p2,rtemp(i)%p3,rtemp(i)%p4,rtemp(i)%p5,rtemp(i)%exothermicity,rtemp(i)%exothermicity_known, rtemp(i)%gamma !(des_reactive*P)/(1+des_reactive*P)
    ENDDO
    CLOSE(1001)


    ! Now resize the reactions array
    DEALLOCATE( r )
    ALLOCATE( r(SIZE(rtemp)) )
    r = rtemp
    DEALLOCATE( rtemp )

    ! Update the number of reactions
    nreactions = nreactions + Nlines

    ! Now read the class 2 suprathermal reactions network
    io     = 0
    Nlines = 0
    OPEN(3,FILE='class_2_suprathermal.dat',STATUS='OLD',IOSTAT=io)
    DO
      READ(3,*,IOSTAT=io)
      IF (io .NE. 0) EXIT
      Nlines = Nlines + 1
    ENDDO
    CLOSE(3)

    ! Allocate temp reactions object to hold the new suprathermal
    ALLOCATE( rtemp(nreactions + Nlines) )
    rtemp(1:SIZE(r)) = r

    ! Add the new suprathermal processes to the temp reactions object
    OPEN(3,FILE='class_2_suprathermal.dat',STATUS='OLD',IOSTAT=io)
    DO ii=nreactions+1,nreactions+Nlines
      READ(3,1000)rtemp(ii)%idx, rtemp(ii)%r1, rtemp(ii)%r2, rtemp(ii)%p1, rtemp(ii)%p2, rtemp(ii)%p3, rtemp(ii)%p4, rtemp(ii)%p5, rtemp(ii)%alpha, rtemp(ii)%beta, rtemp(ii)%gamma
      rtemp(ii)%idx = ii
      rtemp(ii)%ir1 = species_idx(rtemp(ii)%r1)
      rtemp(ii)%ir2 = species_idx(rtemp(ii)%r2)
      rtemp(ii)%ip1 = species_idx(rtemp(ii)%p1)
      rtemp(ii)%ip2 = species_idx(rtemp(ii)%p2)
      rtemp(ii)%ip3 = species_idx(rtemp(ii)%p3)
      rtemp(ii)%ip4 = species_idx(rtemp(ii)%p4)
      rtemp(ii)%ip5 = species_idx(rtemp(ii)%p5)

      rtemp(ii)%exothermicity_known = 0
      rtemp(ii)%exothermicity = 0.0e0
      IF (rtemp(ii)%ir1==-1 .OR. rtemp(ii)%ir2==-1 .OR. rtemp(ii)%ip1==-1 .OR. rtemp(ii)%ip2==-1 .OR. rtemp(ii)%ip3==-1 .OR. rtemp(ii)%ip4==-1 .OR. rtemp(ii)%ip5==-1) THEN
          PRINT '(a31, i5, 1x, 7a10)','Reaction with unknown species: ', rtemp(ii)%idx, rtemp(ii)%r1, rtemp(ii)%r2, rtemp(ii)%p1, rtemp(ii)%p2, rtemp(ii)%p3, rtemp(ii)%p4, rtemp(ii)%p5
          STOP
      ENDIF


      rtemp(ii)%rtype = get_rtype(rtemp(ii)%r1,rtemp(ii)%r2)
    END DO
    CLOSE(3)

    ! Now resize the reactions array
    DEALLOCATE( r )
    ALLOCATE( r(SIZE(rtemp)) )
    r = rtemp
    DEALLOCATE( rtemp )


    ! Print out radiolysis reactions
    OPEN(1001,FILE="class2_reactions.out",STATUS='REPLACE')
    DO i=nreactions+1,nreactions + Nlines
        write(1001,1001)r(i)%idx,r(i)%rtype,r(i)%r1,r(i)%r2,r(i)%p1,r(i)%p2,r(i)%p3,r(i)%p4,r(i)%p5,r(i)%exothermicity,r(i)%exothermicity_known, r(i)%gamma !(des_reactive*P)/(1+des_reactive*P)
    ENDDO
    CLOSE(1001)

    ! Update the number of reactions
    nreactions = nreactions + Nlines

    ! Add quenching reactions for bulk
    ! 1) Allocate temp reactions object to hold the quenching reactions
    !    NB: There should be nspecies - first_suprathermal_species - 1 new
    !    quenching reactions
    ALLOCATE( rtemp(nreactions + (nspecies - first_suprathermal_species + 1) ) )
    rtemp(1:SIZE(r)) = r

    ! Add the new quenching reactions to the temp reactions object
    OPEN(3,FILE='quenching_bulk.out',STATUS='REPLACE',IOSTAT=io)
    jj = first_suprathermal_species
    PRINT *, "First suprathermal_species =",jj
    PRINT *, "nspecies - first_suprathermal_species + 1 = ", nspecies - first_suprathermal_species + 1
    PRINT *, "Species before first suprathermal species =",s(jj-1)%name
    DO ii=nreactions+1,nreactions+(nspecies - first_suprathermal_species + 1)
      groundstate = s(jj)%name(1:LEN_TRIM(s(jj)%name)-1)
      rtemp(ii)%r1  = s(jj)%name
      rtemp(ii)%r2  = "QUENCH"
      rtemp(ii)%p1  = s(jj)%name(1:LEN_TRIM(s(jj)%name)-1)
      rtemp(ii)%p2  = " "
      rtemp(ii)%p3  = " "
      rtemp(ii)%p4  = " "
      rtemp(ii)%p5  = " "
      rtemp(ii)%idx = ii
      rtemp(ii)%ir1 = species_idx(s(jj)%name)
      rtemp(ii)%ir2 = 0
      rtemp(ii)%ip1 = species_idx(groundstate)
      rtemp(ii)%ip2 = 0
      rtemp(ii)%ip3 = 0
      rtemp(ii)%ip4 = 0
      rtemp(ii)%alpha = 1.00e0
      rtemp(ii)%beta = 1.00e0
      rtemp(ii)%gamma = 1.00e0
      rtemp(ii)%rtype = get_rtype(rtemp(ii)%r1,rtemp(ii)%r2)
      rtemp(ii)%exothermicity = 0.00e0
      rtemp(ii)%exothermicity_known = 0
      WRITE(3,1000)rtemp(ii)%idx, rtemp(ii)%r1, rtemp(ii)%r2, rtemp(ii)%p1, rtemp(ii)%p2, rtemp(ii)%p3, rtemp(ii)%p4, rtemp(ii)%p5, rtemp(ii)%alpha, rtemp(ii)%beta, rtemp(ii)%gamma
      PRINT *, s(jj)%name," -> ",s(jj)%name(1:LEN_TRIM(s(jj)%name)-1)
      PRINT *, "jj = ",jj, " nspecies= ",nspecies
      PRINT *, "************************"
      jj = jj + 1
    END DO
    CLOSE(3)

    ! Now resize the reactions array
    DEALLOCATE( r )
    ALLOCATE( r(SIZE(rtemp)) )
    r = rtemp
    DEALLOCATE( rtemp )

    ! Update the number of reactions
    nreactions = nreactions+ (nspecies - first_suprathermal_species + 1)

    PRINT *, s(first_surface_suprathermal)%name
    PRINT *, first_surface_suprathermal
    PRINT *, s(first_suprathermal_species)%name
    PRINT *, s(first_suprathermal_species-1)%name


    ! Add quenching reactions for surface
    ! 1) Allocate temp reactions object to hold the quenching reactions
    !    NB: There should be nspecies - first_suprathermal_species - 1 new
    !    quenching reactions
    ALLOCATE( rtemp(nreactions + (nspecies - first_suprathermal_species + 1) ) )
    rtemp(1:SIZE(r)) = r

    ! Add the new quenching reactions to the temp reactions object
    OPEN(3,FILE='quenching_surface.out',STATUS='replace',IOSTAT=io)
    jj = first_suprathermal_species
    PRINT *, "First suprathermal_species =",jj
    PRINT *, "Species before first suprathermal species =",s(jj-1)%name
    DO ii=nreactions+1,nreactions+(nspecies - first_suprathermal_species + 1)
      groundstate = s(jj)%name(1:LEN_TRIM(s(jj)%name)-1)
      rtemp(ii)%r1  = "g" // s(jj)%name(2:)
      rtemp(ii)%r2  = "QUENCH"
      rtemp(ii)%p1  = "g" // s(jj)%name(2:LEN_TRIM(s(jj)%name)-1)
      rtemp(ii)%p2  = " "
      rtemp(ii)%p3  = " "
      rtemp(ii)%p4  = " "
      rtemp(ii)%p5  = " "
      rtemp(ii)%idx = ii
      rtemp(ii)%ir1 = species_idx(s(jj)%name)
      rtemp(ii)%ir2 = 0
      rtemp(ii)%ip1 = species_idx(groundstate)
      rtemp(ii)%ip2 = 0
      rtemp(ii)%ip3 = 0
      rtemp(ii)%ip4 = 0
      rtemp(ii)%alpha = 1.00e0
      rtemp(ii)%beta = 1.00e0
      rtemp(ii)%gamma = 1.00e0
      rtemp(ii)%rtype = get_rtype(rtemp(ii)%r1,rtemp(ii)%r2)
      rtemp(ii)%exothermicity = 0.00e0
      rtemp(ii)%exothermicity_known = 0
      WRITE(3,1000)rtemp(ii)%idx, rtemp(ii)%r1, rtemp(ii)%r2, rtemp(ii)%p1, rtemp(ii)%p2, rtemp(ii)%p3, rtemp(ii)%p4, rtemp(ii)%p5, rtemp(ii)%alpha, rtemp(ii)%beta, rtemp(ii)%gamma
      PRINT *, "g" // s(jj)%name(2:)," -> g",s(jj)%name(2:LEN_TRIM(s(jj)%name)-1)
      PRINT *, "************************"
      jj = jj + 1
    END DO
    CLOSE(3)

    ! Now resize the reactions array
    DEALLOCATE( r )
    ALLOCATE( r(SIZE(rtemp)) )
    r = rtemp
    DEALLOCATE( rtemp )

    ! Update the number of reactions
    nreactions = nreactions+ (nspecies - first_suprathermal_species + 1)

    DO i=1,nreactions
      PRINT *, r(i)
    END DO

    !***************************************************************************
    ! BEGIN ADD NEW PHOTOCHEMISTRY
    !***************************************************************************

    ! Now read file with radiolysis reactions
    ! and count the number of processes
    io     = 0
    Nlines = 0
    OPEN(3,FILE='photo_processes_2.dat',STATUS='OLD',IOSTAT=io)
    DO
      READ(3,*,IOSTAT=io)
      IF (io .NE. 0) EXIT
      Nlines = Nlines + 1
    ENDDO
    CLOSE(3)

    ! Allocate temp reactions object to hold the radiolysis processes
    ALLOCATE( rtemp(nreactions + Nlines) )
    rtemp(1:SIZE(r)) = r

    ! Add the radiolysis processes to the temp reactions object
    OPEN(3,FILE='photo_processes_2.dat',STATUS='OLD',IOSTAT=io)
    OPEN(1001,FILE="photochemistry_reactions.out",STATUS='REPLACE')
    DO ii=nreactions+1,nreactions + Nlines
      READ(3,1002)rtemp(ii)%idx, rtemp(ii)%r1, rtemp(ii)%r2, rtemp(ii)%p1, rtemp(ii)%p2, rtemp(ii)%p3, rtemp(ii)%p4, rtemp(ii)%p5, rtemp(ii)%alpha, rtemp(ii)%beta, rtemp(ii)%gamma
      WRITE(1001,1002)rtemp(ii)%idx, rtemp(ii)%r1, rtemp(ii)%r2, rtemp(ii)%p1, rtemp(ii)%p2, rtemp(ii)%p3, rtemp(ii)%p4, rtemp(ii)%p5, rtemp(ii)%alpha, rtemp(ii)%beta, rtemp(ii)%gamma
      rtemp(ii)%idx = ii
      rtemp(ii)%ir1 = species_idx(rtemp(ii)%r1)
      rtemp(ii)%ir2 = species_idx(rtemp(ii)%r2)
      rtemp(ii)%ip1 = species_idx(rtemp(ii)%p1)
      rtemp(ii)%ip2 = species_idx(rtemp(ii)%p2)
      rtemp(ii)%ip3 = species_idx(rtemp(ii)%p3)
      rtemp(ii)%ip4 = species_idx(rtemp(ii)%p4)
      rtemp(ii)%ip5 = species_idx(rtemp(ii)%p5)
      IF (rtemp(ii)%ir1==-1 .OR. rtemp(ii)%ir2==-1 .OR. rtemp(ii)%ip1==-1 .OR. rtemp(ii)%ip2==-1 .OR. rtemp(ii)%ip3==-1 .OR. rtemp(ii)%ip4==-1 .OR. rtemp(ii)%ip5==-1) THEN
          PRINT '(a31, i5, 1x, 7a10)','Reaction with unknown species: ', rtemp(ii)%idx, rtemp(ii)%r1, rtemp(ii)%r2, rtemp(ii)%p1, rtemp(ii)%p2, rtemp(ii)%p3, rtemp(ii)%p4, rtemp(ii)%p5
          STOP
      ENDIF

      rtemp(ii)%rtype = get_rtype(rtemp(ii)%r1,rtemp(ii)%r2)
      rtemp(ii)%exothermicity_known = 0.0
      rtemp(ii)%exothermicity = 0.0
    END DO
    CLOSE(3)
    CLOSE(1001)


    ! Now resize the reactions array
    DEALLOCATE( r )
    ALLOCATE( r(SIZE(rtemp)) )
    r = rtemp
    DEALLOCATE( rtemp )

    ! Update the number of reactions
    nreactions = nreactions + Nlines
    !***************************************************************************
    ! END ADD NEW PHOTOCHEMISTRY
    !***************************************************************************

    ! Now loop through all reactions and make sure they are correct
    DO i=1,nreactions
      IF ( r(i)%ir1 .NE. 0 ) THEN
        ! Sum prod atoms
        prodatoms = s(r(i)%ip1)%natoms
        IF ( r(i)%ip2 .NE. 0 ) prodatoms = prodatoms + s(r(i)%ip2)%natoms
        IF ( r(i)%ip3 .NE. 0 ) prodatoms = prodatoms + s(r(i)%ip3)%natoms
        IF ( r(i)%ip4 .NE. 0 ) prodatoms = prodatoms + s(r(i)%ip4)%natoms
        reactatoms = s(r(i)%ir1)%natoms

        IF (ANY(r(i)%r2 .EQ. (/"QUENCH","CRPHOT","PHOTON","FREEZE","DESORB","IONRAD","G-    ","G0    ","CR    ","CRP   ","PHOION","PHOEXC"/)) .EQV. .FALSE.) THEN
          reactatoms = reactatoms + s(r(i)%ir2)%natoms
        ENDIF

        PRINT *, "-----------------------------------------"
        PRINT *, r(i)%r1," + ",r(i)%r2," -> ",r(i)%p1," + ",r(i)%p2
        PRINT *, "Prodatoms=",prodatoms
        PRINT *, "Reactatoms=",reactatoms
        PRINT *, "ir1_natoms = ",s(r(i)%ir1)%natoms
        PRINT *, "ip1_natoms = ",s(r(i)%ip1)%natoms
        PRINT *, "alpha=",r(i)%alpha," beta=",r(i)%beta," gamma=",r(i)%gamma
        PRINT *, "rate=",r(i)%rate, " exothermicity =",r(i)%exothermicity
        PRINT *, r(i)%ir1,r(i)%ir2,r(i)%ip1,r(i)%ip2,r(i)%ip3,r(i)%ip4
        IF (prodatoms .NE. reactatoms) THEN
          IF  ( r(i)%r1(2:2) .NE. "e" ) THEN
            PRINT *, "Atoms not equal"
            CALL EXIT()
          ENDIF
        ENDIF
        IF ( ISNAN(r(i)%alpha) ) THEN
          PRINT *, "Alpha = NaN"
        ENDIF

        IF ( ISNAN(r(i)%beta) ) THEN
          PRINT *, "Beta = NaN"
        ENDIF

        IF ( ISNAN(r(i)%gamma) ) THEN
          PRINT *, "Gamma = NaN"
        ENDIF

        IF ( ISNAN(r(i)%rate) ) THEN
          PRINT *, "Rate = NaN"
        ENDIF

        IF ( ISNAN(r(i)%exothermicity) ) THEN
          PRINT *, "Exothermicity = NaN"
        ENDIF


      ENDIF
    ENDDO
END IF



OPEN (1,FILE='species.out', STATUS='REPLACE')
write(1,*)'Species, Weight, Number of Atoms, Index, Gas twin index, Desorption energy'
DO i = 1, nspecies
  write(1,'(a10,1x,1pe12.4,1x,i3,1x,2i4,1pe12.4)')s(i)%name, s(i)%weight, s(i)%natoms, s(i)%idx, s(i)%gas_idx, s(i)%edes
ENDDO
CLOSE (1)

RETURN
100 PRINT*, 'File ',chem_file(1:LEN_TRIM(chem_file)),' not found!'
STOP
1000 FORMAT(1X,I4,1X,2(A10),10X,5(A10),E8.2,1X,F5.2,1X,F8.1)
1001 FORMAT(1X,I4,1X,I2,1X,2(A10),10X,5(A10),1pE14.5,1X,I1,1X,1pE14.5)
1002 FORMAT(1X,I4,1X,2(A10),10X,5(A10),ES8.2,2X,ES8.2,2X,ES8.2)
END SUBROUTINE read_rate06database

SUBROUTINE read_enthalpias
IMPLICIT NONE
CHARACTER*10 :: local_species
REAL*8       :: local_enthalpia
INTEGER      :: local_enthalpia_known, i

OPEN(2, FILE='enthalpias.txt', STATUS='OLD', ERR=200)

READ(2,*)
DO WHILE (.not. EOF(unit=2) )
    READ(2, '(a10,1x,e8.2,2x,i1)')local_species, local_enthalpia, local_enthalpia_known
    DO i = 1, nspecies
        IF (s(i)%name=='g'//local_species .OR. s(i)%name==local_species) THEN
            s(i)%enthalpia = local_enthalpia * 1.0d3/8.31
            s(i)%enthalpia_known = local_enthalpia_known
        ENDIF
    ENDDO
ENDDO
CLOSE(2)

RETURN
200 PRINT*, 'File enthalpias.txt not found!'
STOP
END SUBROUTINE read_enthalpias

SUBROUTINE get_reaction_thermodynamics(ii)
IMPLICIT NONE

INTEGER :: ii

  !Exothermicity of reaction
  r(ii)%exothermicity_known = s(r(ii)%ir1)%enthalpia_known
  IF (r(ii)%ir2/=0) r(ii)%exothermicity_known = r(ii)%exothermicity_known*s(r(ii)%ir2)%enthalpia_known
  IF (r(ii)%ip1/=0) r(ii)%exothermicity_known = r(ii)%exothermicity_known*s(r(ii)%ip1)%enthalpia_known
  IF (r(ii)%ip2/=0) r(ii)%exothermicity_known = r(ii)%exothermicity_known*s(r(ii)%ip2)%enthalpia_known
  IF (r(ii)%ip3/=0) r(ii)%exothermicity_known = r(ii)%exothermicity_known*s(r(ii)%ip3)%enthalpia_known
  IF (r(ii)%ip4/=0) r(ii)%exothermicity_known = r(ii)%exothermicity_known*s(r(ii)%ip4)%enthalpia_known
  IF (r(ii)%ip5/=0) r(ii)%exothermicity_known = r(ii)%exothermicity_known*s(r(ii)%ip5)%enthalpia_known

    r(ii)%exothermicity = 0.0d0
  IF (r(ii)%exothermicity_known==1) THEN
        r(ii)%exothermicity = r(ii)%exothermicity - s(r(ii)%ir1)%enthalpia
      IF (r(ii)%ir2/=0) r(ii)%exothermicity = r(ii)%exothermicity - s(r(ii)%ir2)%enthalpia
      IF (r(ii)%ip1/=0) r(ii)%exothermicity = r(ii)%exothermicity + s(r(ii)%ip1)%enthalpia
      IF (r(ii)%ip2/=0) r(ii)%exothermicity = r(ii)%exothermicity + s(r(ii)%ip2)%enthalpia
      IF (r(ii)%ip3/=0) r(ii)%exothermicity = r(ii)%exothermicity + s(r(ii)%ip3)%enthalpia
      IF (r(ii)%ip4/=0) r(ii)%exothermicity = r(ii)%exothermicity + s(r(ii)%ip4)%enthalpia
      IF (r(ii)%ip5/=0) r(ii)%exothermicity = r(ii)%exothermicity + s(r(ii)%ip5)%enthalpia
    ENDIF

  r(ii)%exothermicity = -r(ii)%exothermicity

RETURN
STOP
END SUBROUTINE get_reaction_thermodynamics

SUBROUTINE get_rd_efficiency(ii)
IMPLICIT NONE

INTEGER :: ii
REAL*8  :: P

  !Setting efficiency of reactive desorption using des_reactive_type
    SELECT CASE (des_reactive_type)
    CASE (0) !Garrod_ea07
      IF (r(ii)%rtype == 13 .AND. r(ii)%p1(1:1) == 'g') THEN
          IF (r(ii)%ip2/=0) THEN
              r(ii)%alpha = 1.0d0*r(ii)%alpha
          ELSE
              IF (r(ii)%exothermicity_known==1) THEN
                  P = (1.d0-s(r(ii)%ip1)%edes/r(ii)%exothermicity)**(dmax1(2.d0,3.d0*s(r(ii)%ip1)%natoms-5.d0)-1.d0)
                  r(ii)%alpha =(1.d0 - (des_reactive*P)/(1+des_reactive*P))*r(ii)%alpha
              ELSE
                  r(ii)%alpha = (1.d0 - des_reactive)*r(ii)%alpha
              ENDIF
          ENDIF
      ENDIF
      IF (r(ii)%rtype == 13 .AND. r(ii)%p1(1:1) /= 'g') THEN
          IF (r(ii)%ip2/=0) THEN
              r(ii)%alpha = 0.d0
          ELSE
              IF (r(ii)%exothermicity_known==1) THEN
                  P = (1.d0-s(r(ii)%ip1)%edes/r(ii)%exothermicity)**(dmax1(2.d0,3.d0*s(r(ii)%ip1)%natoms-5.d0)-1.d0)
                  r(ii)%alpha =((des_reactive*P)/(1+des_reactive*P))*r(ii)%alpha
              ELSE
                  r(ii)%alpha = des_reactive*r(ii)%alpha
              ENDIF
          ENDIF
      ENDIF
  CASE (1) !Vasyunin&Herbst13
      IF (r(ii)%rtype == 13 .AND. r(ii)%p1(1:1) == 'g') r(ii)%alpha = (1.d0 - des_reactive)*r(ii)%alpha
      IF (r(ii)%rtype == 13 .AND. r(ii)%p1(1:1) /= 'g') r(ii)%alpha = des_reactive*r(ii)%alpha
  CASE (2) !Minissale&Dulieu
      !!!des_reactive = 0.0d0
      IF (r(ii)%rtype == 13 .AND. r(ii)%p1(1:1) == 'g') r(ii)%alpha = (1.d0 - des_reactive)*r(ii)%alpha
      IF (r(ii)%rtype == 13 .AND. r(ii)%p1(1:1) /= 'g') r(ii)%alpha = des_reactive*r(ii)%alpha
      !Desorption efficiency is individual for each species, not reaction! As such, it is calculated in mod_run_dvode.f90
  CASE DEFAULT
      PRINT*, 'Unknown treatment of reactive desorption: ',des_reactive_type
      STOP
  END SELECT

RETURN
STOP
END SUBROUTINE get_rd_efficiency

!This function returns the reactions type by its reactants
!
!Reaction types:
!
! 1.  Two-body molecule(ion) - molecule(ion) gas-phase reaction
! 2.  Cosmic ray ionization reaction (CRP)
! 3.  Photoionization reaction (PHOTON)
! 4.  Photoionization, H2 self-shielded
! 5.  Photoionization, CO self-shielded
! 6.  Cosmic ray-induced photoreaction (CRPHOT)
! 7.  Ions plus negatively charged grains
! 8.  Ions plus neutral grains
! 9.  Collisions of the electrons with the neutral grains
! 10. Collisions of the electrons with the positively charged grains
! 11. Accretion (FREEZE)
! 12. Desorption (DESORB)
! 13. Surface two-body reaction
! 14. Bulk two-body reaction
! 15. Suprathermal surface reaction
! 16. Suprathermal bulk reaction
! 17. Solid-phase radiolysis (IONRAD)
! 18. Quenching of suprathermal species (QUENCH)
! 19. Photoionization
! 20. Photoexcitation
! 21. Surface ion-neutral
! 22. Bulk ion-neutral
INTEGER FUNCTION get_rtype(r1,r2)
IMPLICIT NONE
CHARACTER*10 r1, r2

get_rtype = 1

IF (r2(1:LEN_TRIM(r2)) == 'CRP') get_rtype = 2
IF (r2(1:LEN_TRIM(r2)) == 'PHOTON') get_rtype = 3
IF (r1(1:LEN_TRIM(r1)) == 'H2' .AND. r2(1:LEN_TRIM(r2)) == 'PHOTON') get_rtype = 4
IF (r1(1:LEN_TRIM(r1)) == 'CO' .AND. r2(1:LEN_TRIM(r2)) == 'PHOTON') get_rtype = 5
IF (r2(1:LEN_TRIM(r2)) == 'CRPHOT') get_rtype = 6
IF (r2(1:LEN_TRIM(r2)) == 'G-') get_rtype = 7
IF (r2(1:LEN_TRIM(r2)) == 'G0') get_rtype = 8
IF (r1(1:LEN_TRIM(r1)) == 'G0') get_rtype = 9
IF (r1(1:LEN_TRIM(r1)) == 'G+') get_rtype = 10
IF (r2(1:LEN_TRIM(r2)) == 'FREEZE') get_rtype = 11
IF (r2(1:LEN_TRIM(r2)) == 'DESORB') get_rtype = 12
IF (r1(1:1) == 'g' .AND. r2(1:1) == 'g') get_rtype = 13
IF (r1(1:1) == 'b' .AND. r2(1:1) == 'b') get_rtype = 14

IF (r1(LEN_TRIM(r1):LEN_TRIM(r1)) .EQ. '+' .OR. r1(LEN_TRIM(r1):LEN_TRIM(r1)) .EQ. '-' ) THEN 
  IF ( r1(1:1) .EQ. 'g' ) THEN
    get_rtype = 21
  ELSE IF (r1(1:1).EQ.'b' ) THEN
    get_rtype = 22
  ENDIF
ENDIF

IF (r2(LEN_TRIM(r2):LEN_TRIM(r2)) .EQ. '+' .OR. r2(LEN_TRIM(r2):LEN_TRIM(r2)) .EQ. '-' ) THEN
  IF ( r1(1:1) .EQ. 'g' ) THEN
    get_rtype = 21
  ELSE IF (r1(1:1).EQ.'b' ) THEN
    get_rtype = 22
  ENDIF
ENDIF

IF (((r1(LEN_TRIM(r1):LEN_TRIM(r1)) .EQ. '*' .OR. r2(LEN_TRIM(r2):LEN_TRIM(r2)) .EQ. '*')) &
  .AND. ((r1(1:1) .EQ. 'g') .OR. (r2(1:1) .EQ. 'g'))) get_rtype = 15
IF (((r1(LEN_TRIM(r1):LEN_TRIM(r1)) .EQ. '*' .OR. r2(LEN_TRIM(r2):LEN_TRIM(r2)) .EQ. '*')) &
  .AND. ((r1(1:1) .EQ. 'b') .OR. (r2(1:1) .EQ. 'b'))) get_rtype = 16
IF (r2(1:LEN_TRIM(r2)) == 'IONRAD') get_rtype = 17
IF (r2(1:LEN_TRIM(r2)) == 'QUENCH') THEN
  PRINT *, r1," + ",r2
  get_rtype = 18
ENDIF
IF (r2(1:LEN_TRIM(r2)) == 'PHOION') get_rtype = 19
IF (r2(1:LEN_TRIM(r2)) == 'PHOEXC') get_rtype = 20

END FUNCTION get_rtype

END
