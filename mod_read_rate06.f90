MODULE read_rate06
  USE global_variables
  USE global_functions
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_rate06database, read_enthalpias, get_reaction_thermodynamics, &
            get_rd_efficiency, get_rtype

  INTEGER, PARAMETER :: MAX_MONOLAYERS_DEFAULT = 10000
  INTEGER, PARAMETER :: MIN_UNIT = 10, MAX_UNIT = 999
  INTEGER, PARAMETER :: MAX_FILE_PATH_LEN = 255

CONTAINS

SUBROUTINE read_rate06database
  IMPLICIT NONE
  INTEGER :: i, ii, jj, idx, s_r_counter, Nsup_g, Nlines, io_stat, file_unit
  INTEGER :: prodatoms, reactatoms, alloc_stat, num_fields
  INTEGER :: first_suprathermal_react, first_suprathermal_species, max_monolayers
  INTEGER :: output_unit
  INTEGER :: j, k  ! Loop variables
  REAL(KIND=wp) :: a, b, c, apriori_nml
  CHARACTER(LEN=10) :: groundstate, s_name, r1, r2, p1, p2, p3, p4, p5
  CHARACTER(LEN=MAX_FILE_PATH_LEN) :: filename
  TYPE(reaction), DIMENSION(:), ALLOCATABLE :: rtemp
  LOGICAL :: file_exists, is_suprathermal_on
  CHARACTER(LEN=120) :: line
  CHARACTER(LEN=10) :: temp_fields(8)  ! Temporary array for all fields after reactants

  ! Initialization
  nreactions = 0
  nspecies = 0
  n_surf_spec = 0
  n_surf_react = 0
  first_surf_spec = 0
  first_surfreact = 0
  first_bulkreact = 0
  first_suprathermal_species = 0
  Nsup_g = 0
  first_suprathermal_react = 0
  file_unit = MIN_UNIT
  output_unit = MIN_UNIT

  ! Set max_monolayers
  IF (MODEL_EXPERIMENT == 1) THEN
    apriori_nml = MIN(1.0e6_wp, MAX(100.0_wp, 10.0_wp * (ICE_THICK / 5.0e-8_wp)))
  ELSE
    apriori_nml = 1000.0_wp
  END IF
  max_monolayers = MIN(MAX_MONOLAYERS_DEFAULT, INT(apriori_nml))

  ! File check with debug
  WRITE(*,'(A,A)') 'chem_file = ', TRIM(chem_file)
  WRITE(*,'(A,I0)') 'Length of chem_file = ', LEN_TRIM(chem_file)
  INQUIRE(FILE=TRIM(chem_file), EXIST=file_exists)
  IF (.NOT. file_exists) THEN
    WRITE(*,'(A,A,A)') 'ERROR: File "', TRIM(chem_file), '" not found!'
    CALL cleanup_and_exit()
  END IF

  ! Find unit and open file
  WRITE(*,'(A)') 'Calling find_free_unit'
  CALL find_free_unit(file_unit, MIN_UNIT, MAX_UNIT)
  WRITE(*,'(A,I0)') 'file_unit = ', file_unit
  IF (file_unit == -1) THEN
    WRITE(*,*) 'ERROR: No available file units'
    CALL cleanup_and_exit()
  END IF
  WRITE(*,'(A)') 'Opening file for pre-read and main read'
  OPEN(UNIT=file_unit, FILE=TRIM(chem_file), STATUS='OLD', IOSTAT=io_stat)
  IF (io_stat /= 0) THEN
    WRITE(*,'(A,I0)') 'ERROR: Failed to open file, IOSTAT=', io_stat
    CALL cleanup_and_exit()
  END IF

  ! Pre-read
  WRITE(*,'(A)') 'Calling pre_read_database'
  CALL pre_read_database(file_unit)
  WRITE(*,'(A)') 'Finished pre_read_database'
  REWIND(file_unit)
  WRITE(*,'(A)') 'File rewound for main read'

  ! Allocate arrays with full anticipated size
  is_suprathermal_on = (suprathermal == 1)
  WRITE(*,'(A,I0)') 'Allocating arrays with nspecies=', nspecies
  ! Anticipate max species: initial (22) + suprathermal surface (10) + bulk (10) + suprathermal bulk (10)
  IF (is_suprathermal_on) THEN
    CALL allocate_arrays_with_size(52, max_monolayers)  ! 52 = 22 + 10 + 10 + 10
  ELSE
    CALL allocate_arrays_with_size(42, max_monolayers)  ! 42 = 22 + 10 + 10
  END IF

  ! Main reading phase - species first
  WRITE(*,'(A)') 'Starting main read'
  READ(file_unit, *, IOSTAT=io_stat) nspecies
  IF (io_stat /= 0) CALL handle_io_error('reading nspecies')
  DO i = 1, nspecies
    READ(file_unit, '(A10)', IOSTAT=io_stat) s(i)%name
    IF (io_stat /= 0) CALL handle_io_error('reading species name', i)
    s(i)%idx = i
    s(i)%gas_idx = i
    s(i)%weight = aweight(s(i)%name)
    s(i)%natoms = numatoms(s(i)%name)
    IF (first_surf_spec == 0 .AND. s(i)%name(1:1) == 'g') first_surf_spec = i
  END DO

  ! Link gas counterparts
  CALL link_gas_counterparts()

  ! Add suprathermal surface species AFTER main species read
  IF (is_suprathermal_on) THEN
    WRITE(*,'(A)') 'Adding suprathermal surface species after main species read'
    CALL add_suprathermal_surface_species()
    ! Debug: List all species after adding suprathermal surface species
    WRITE(*,'(A,I0)') 'Total species after suprathermal surface addition: ', nspecies
    DO i = 1, nspecies
      WRITE(*,'(A,I0,A,A10)') 'Species ', i, ': ', TRIM(s(i)%name)
    END DO
  END IF

  ! Read enthalpies
  CALL read_enthalpias()

  ! Read reactions
  READ(file_unit, *, IOSTAT=io_stat) nreactions
  IF (io_stat /= 0) CALL handle_io_error('reading nreactions')
  CALL allocate_reaction_arrays(is_suprathermal_on)
  ii = 1
  DO i = 1, nreactions
    READ(file_unit, '(A)', IOSTAT=io_stat) line
    IF (io_stat /= 0) THEN
      WRITE(*,'(A,I0,A,I0)') 'Failed to read reaction line ', i, ' with IOSTAT=', io_stat
      CALL handle_io_error('reading reaction raw', i)
    END IF
    WRITE(*,'(A,I0,A,A)') 'Reaction line ', i, ': ', TRIM(line)
    ! Initial parse for core fields
    READ(line, *, IOSTAT=io_stat) r(ii)%idx, r(ii)%r1, r(ii)%r2
    IF (io_stat /= 0) THEN
      WRITE(*,'(A,I0,A,I0)') 'Failed to parse reaction core ', i, ' with IOSTAT=', io_stat
      CALL handle_io_error('reading reaction core', i)
    END IF
    ! Read all fields into temp array
    temp_fields = ' '
    READ(line, *, IOSTAT=io_stat) r(ii)%idx, r(ii)%r1, r(ii)%r2, &
      temp_fields(1), temp_fields(2), temp_fields(3), temp_fields(4), temp_fields(5), &
      temp_fields(6), temp_fields(7), temp_fields(8)
    IF (io_stat /= 0 .AND. io_stat /= -1) THEN
      WRITE(*,'(A,I0,A,I0)') 'Failed to parse reaction fields ', i, ' with IOSTAT=', io_stat
      CALL handle_io_error('reading reaction fields', i)
    END IF
    ! Determine number of fields after reactants
    num_fields = 0
    DO j = 1, 8
      IF (LEN_TRIM(temp_fields(j)) > 0) THEN
        num_fields = num_fields + 1
      ELSE
        EXIT
      END IF
    END DO
    ! Assign products and reals based on field count
    IF (num_fields >= 4 .AND. num_fields <= 8) THEN  ! 1 to 5 products + 3 reals
      SELECT CASE (num_fields)
      CASE (4)  ! 1 product
        r(ii)%p1 = temp_fields(1)
        r(ii)%p2 = ' '
        r(ii)%p3 = ' '
        r(ii)%p4 = ' '
        r(ii)%p5 = ' '
        r(ii)%alpha = REALVALUE(temp_fields(2))
        r(ii)%beta = REALVALUE(temp_fields(3))
        r(ii)%gamma = REALVALUE(temp_fields(4))
      CASE (5)  ! 2 products
        r(ii)%p1 = temp_fields(1)
        r(ii)%p2 = temp_fields(2)
        r(ii)%p3 = ' '
        r(ii)%p4 = ' '
        r(ii)%p5 = ' '
        r(ii)%alpha = REALVALUE(temp_fields(3))
        r(ii)%beta = REALVALUE(temp_fields(4))
        r(ii)%gamma = REALVALUE(temp_fields(5))
      CASE (6)  ! 3 products
        r(ii)%p1 = temp_fields(1)
        r(ii)%p2 = temp_fields(2)
        r(ii)%p3 = temp_fields(3)
        r(ii)%p4 = ' '
        r(ii)%p5 = ' '
        r(ii)%alpha = REALVALUE(temp_fields(4))
        r(ii)%beta = REALVALUE(temp_fields(5))
        r(ii)%gamma = REALVALUE(temp_fields(6))
      CASE (7)  ! 4 products
        r(ii)%p1 = temp_fields(1)
        r(ii)%p2 = temp_fields(2)
        r(ii)%p3 = temp_fields(3)
        r(ii)%p4 = temp_fields(4)
        r(ii)%p5 = ' '
        r(ii)%alpha = REALVALUE(temp_fields(5))
        r(ii)%beta = REALVALUE(temp_fields(6))
        r(ii)%gamma = REALVALUE(temp_fields(7))
      CASE (8)  ! 5 products
        r(ii)%p1 = temp_fields(1)
        r(ii)%p2 = temp_fields(2)
        r(ii)%p3 = temp_fields(3)
        r(ii)%p4 = temp_fields(4)
        r(ii)%p5 = temp_fields(5)
        r(ii)%alpha = REALVALUE(temp_fields(6))
        r(ii)%beta = REALVALUE(temp_fields(7))
        r(ii)%gamma = REALVALUE(temp_fields(8))
      END SELECT
    ELSE
      WRITE(*,'(A,I0,A,I0)') 'Invalid field count for reaction ', i, ': ', num_fields
      CALL handle_io_error('invalid field count', i)
    END IF
    WRITE(*,'(A,I0)') 'Processing reaction ', r(ii)%idx
    CALL assign_species_indices(ii)
    CALL get_reaction_thermodynamics(ii)
    r(ii)%rtype = get_rtype(r(ii)%r1, r(ii)%r2)
    IF (r(ii)%rtype == 12 .AND. r(ii)%ir1 > 0) THEN
      s(r(ii)%ir1)%edes = r(ii)%gamma
      IF (s(r(ii)%ir1)%gas_idx > 0) s(s(r(ii)%ir1)%gas_idx)%edes = r(ii)%gamma
    END IF
    IF (first_surfreact == 0 .AND. r(ii)%rtype == 13) first_surfreact = ii
    CALL get_rd_efficiency(ii)
    WRITE(37, 1001) r(ii)%idx, r(ii)%rtype, r(ii)%r1, r(ii)%r2, r(ii)%p1, &
      r(ii)%p2, r(ii)%p3, r(ii)%p4, r(ii)%p5, r(ii)%exothermicity, &
      r(ii)%exothermicity_known, r(ii)%alpha
    ii = ii + 1
  END DO
  CLOSE(file_unit)

  ! Post-processing
  WRITE(*,'(A)') 'Finished reading reactions'
  CALL add_bulk_species(first_suprathermal_species)
  WRITE(*,'(A)') 'Added bulk species'
  ! Debug: List species after bulk addition
  WRITE(*,'(A,I0)') 'Total species after bulk addition: ', nspecies
  DO i = 1, nspecies
    WRITE(*,'(A,I0,A,A10)') 'Species ', i, ': ', TRIM(s(i)%name)
  END DO
  CALL add_bulk_reactions(output_unit)
  WRITE(*,'(A,L1)') 'Added bulk reactions, is_suprathermal_on=', is_suprathermal_on
  IF (is_suprathermal_on) THEN
    CALL add_suprathermal_bulk_and_reactions(Nsup_g, first_suprathermal_react, output_unit)
    ! Debug: List species after suprathermal bulk addition
    WRITE(*,'(A,I0)') 'Total species after suprathermal bulk addition: ', nspecies
    DO i = 1, nspecies
      WRITE(*,'(A,I0,A,A10)') 'Species ', i, ': ', TRIM(s(i)%name)
    END DO
  END IF
  WRITE(*,'(A)') 'Calling write_species_output'
  CALL write_species_output(file_unit)
  WRITE(*,'(A)') 'Finished write_species_output'

  RETURN
1000 FORMAT(1X,I4,1X,2(A10),10X,5(A10),E12.4,1X,F8.2,1X,F8.1)  ! For reference
1001 FORMAT(1X,I4,1X,I2,1X,2(A10),10X,5(A10),1pE14.5,1X,I1,1X,1pE14.5)
END SUBROUTINE read_rate06database

SUBROUTINE allocate_arrays_with_size(total_species, max_monolayers)
  INTEGER, INTENT(IN) :: total_species, max_monolayers
  INTEGER :: i, alloc_stat
  ALLOCATE(s(total_species), STAT=alloc_stat)
  IF (alloc_stat /= 0) CALL handle_alloc_error('species array', total_species)
  IF (MODEL_EXPERIMENT == 0) THEN
    ALLOCATE(abundances_bulk(total_species, MAX_MONOLAYERS_DEFAULT), STAT=alloc_stat)
    ALLOCATE(timesteps_nml(3000), STAT=alloc_stat)
  ELSE
    ALLOCATE(abundances_bulk(total_species, max_monolayers), STAT=alloc_stat)
    ALLOCATE(timesteps_nml(max_monolayers), STAT=alloc_stat)
  END IF
  IF (alloc_stat /= 0) CALL handle_alloc_error('abundances_bulk or timesteps_nml', total_species)
  DO i = 1, total_species
    ALLOCATE(s(i)%abundance_out(timesteps), STAT=alloc_stat)
    IF (alloc_stat /= 0) CALL handle_alloc_error('abundance_out', i)
    s(i)%abundance_out(:) = 0.0_wp
    s(i)%edes = 0.0_wp
    s(i)%racc = 0.0_wp
    s(i)%rdes = 0.0_wp
    s(i)%abundance = 0.0_wp
    s(i)%frac_abundance = 0.0_wp
  END DO
END SUBROUTINE allocate_arrays_with_size

SUBROUTINE assign_species_indices(ii)
  INTEGER, INTENT(IN) :: ii
  INTEGER :: ir1, ir2, ip1, ip2, ip3, ip4, ip5
  ir1 = species_idx_with_autocreate(r(ii)%r1)  ! Change to use the auto-create function
  WRITE(*,'(A,A10,A,I0)') 'Species ', r(ii)%r1, ' index=', ir1
  ir2 = species_idx_with_autocreate(r(ii)%r2)  ! Change to use the auto-create function
  WRITE(*,'(A,A10,A,I0)') 'Species ', r(ii)%r2, ' index=', ir2
  ip1 = species_idx_with_autocreate(r(ii)%p1)  ! Change to use the auto-create function
  WRITE(*,'(A,A10,A,I0)') 'Species ', r(ii)%p1, ' index=', ip1
  ip2 = species_idx_with_autocreate(r(ii)%p2)  ! Change to use the auto-create function
  WRITE(*,'(A,A10,A,I0)') 'Species ', r(ii)%p2, ' index=', ip2
  ip3 = species_idx_with_autocreate(r(ii)%p3)  ! Change to use the auto-create function
  WRITE(*,'(A,A10,A,I0)') 'Species ', r(ii)%p3, ' index=', ip3
  ip4 = species_idx_with_autocreate(r(ii)%p4)  ! Change to use the auto-create function
  WRITE(*,'(A,A10,A,I0)') 'Species ', r(ii)%p4, ' index=', ip4
  ip5 = species_idx_with_autocreate(r(ii)%p5)  ! Change to use the auto-create function
  WRITE(*,'(A,A10,A,I0)') 'Species ', r(ii)%p5, ' index=', ip5
  r(ii)%ir1 = ir1
  r(ii)%ir2 = ir2
  r(ii)%ip1 = ip1
  r(ii)%ip2 = ip2
  r(ii)%ip3 = ip3
  r(ii)%ip4 = ip4
  r(ii)%ip5 = ip5
  IF (ANY([r(ii)%ir1, r(ii)%ir2, r(ii)%ip1, r(ii)%ip2, r(ii)%ip3, r(ii)%ip4, r(ii)%ip5] == -1)) THEN
    WRITE(*,'(A,I0,1X,7A10)') 'ERROR: Unknown species in reaction: ', r(ii)%idx, &
      r(ii)%r1, r(ii)%r2, r(ii)%p1, r(ii)%p2, r(ii)%p3, r(ii)%p4, r(ii)%p5
    CALL cleanup_and_exit()
  END IF
END SUBROUTINE assign_species_indices

  SUBROUTINE pre_read_database(file_unit)
    INTEGER, INTENT(IN) :: file_unit
    INTEGER :: i, io_stat, idx, j, k
    CHARACTER(LEN=10) :: s_name, r1, r2, p1, p2
    CHARACTER(LEN=10) :: temp_fields(5)
    REAL(KIND=wp) :: a, b, c
    CHARACTER(LEN=120) :: line
    INTEGER :: num_fields
    io_stat = 0
    READ(file_unit, *, IOSTAT=io_stat) nspecies
    IF (io_stat /= 0) THEN
      WRITE(*,'(A,I0)') 'Failed to read nspecies with IOSTAT=', io_stat
      CALL cleanup_and_exit()
    END IF
    DO i = 1, nspecies
      READ(file_unit, '(A10)', IOSTAT=io_stat) s_name
      IF (io_stat /= 0) THEN
        WRITE(*,'(A,I0)') 'Species read failed with IOSTAT=', io_stat
        CALL cleanup_and_exit()
      END IF
      IF (s_name(1:1) == 'g') n_surf_spec = n_surf_spec + 1
    END DO
    READ(file_unit, *, IOSTAT=io_stat) nreactions
    IF (io_stat /= 0) THEN
      WRITE(*,'(A,I0)') 'Failed to read nreactions with IOSTAT=', io_stat
      CALL cleanup_and_exit()
    END IF
    DO i = 1, nreactions
      READ(file_unit, '(A)', IOSTAT=io_stat) line
      IF (io_stat /= 0) THEN
        WRITE(*,'(A,I0)') 'Failed to read raw line with IOSTAT=', io_stat
        CALL cleanup_and_exit()
      END IF
      WRITE(*,'(A,I0,A,A)') 'Reading reaction ', i, ': ', TRIM(line)
      READ(line, *, IOSTAT=io_stat) idx, r1, r2
      IF (io_stat /= 0) THEN
        WRITE(*,'(A,I0,A,I0)') 'Initial parse failed for reaction ', i, ' with IOSTAT=', io_stat
        CALL cleanup_and_exit()
      END IF
      num_fields = 0
      READ(line, *, IOSTAT=io_stat) idx, r1, r2, (temp_fields(j), j=1,5)
      IF (io_stat == 0) THEN
        num_fields = 5
      ELSE IF (io_stat == -1) THEN
        DO j = 1, 5
          READ(line, *, IOSTAT=io_stat) idx, r1, r2, (temp_fields(k), k=1,j)
          IF (io_stat == -1) THEN
            num_fields = j - 1
            EXIT
          ELSE IF (io_stat /= 0) THEN
            WRITE(*,'(A,I0,A,I0)') 'Field count failed for reaction ', i, ' with IOSTAT=', io_stat
            CALL cleanup_and_exit()
          END IF
        END DO
      ELSE
        WRITE(*,'(A,I0,A,I0)') 'Temp parse failed for reaction ', i, ' with IOSTAT=', io_stat
        CALL cleanup_and_exit()
      END IF
      IF (num_fields < 3 .OR. num_fields > 5) THEN
        WRITE(*,'(A,I0,A,I0)') 'Invalid field count for reaction ', i, ': ', num_fields
        CALL cleanup_and_exit()
      END IF
      p1 = temp_fields(1)
      IF (num_fields >= 4) THEN
        p2 = temp_fields(2)
        a = REALVALUE(temp_fields(3))
        b = REALVALUE(temp_fields(4))
        IF (num_fields == 5) THEN
          c = REALVALUE(temp_fields(5))
        ELSE
          c = 0.0_wp
        END IF
      ELSE
        p2 = ' '
        a = REALVALUE(temp_fields(2))
        b = REALVALUE(temp_fields(3))
        c = 0.0_wp
      END IF
      IF (r1(1:1) == 'g' .AND. r2 /= 'FREEZE' .AND. r2 /= 'DESORB' .AND. p1(1:1) == 'g') &
        n_surf_react = n_surf_react + 1
    END DO
  END SUBROUTINE pre_read_database

  SUBROUTINE find_free_unit(unit_num, min_unit, max_unit)
    INTEGER, INTENT(OUT) :: unit_num
    INTEGER, INTENT(IN) :: min_unit, max_unit
    LOGICAL :: is_open
    INTEGER :: i
    unit_num = -1
    DO i = min_unit, max_unit
      INQUIRE(UNIT=i, OPENED=is_open)
      IF (.NOT. is_open) THEN
        unit_num = i
        EXIT
      END IF
    END DO
  END SUBROUTINE find_free_unit

  SUBROUTINE allocate_arrays(is_suprathermal_on, max_monolayers)
    LOGICAL, INTENT(IN) :: is_suprathermal_on
    INTEGER, INTENT(IN) :: max_monolayers
    INTEGER :: total_species, i, alloc_stat
    total_species = nspecies + MERGE(3 * n_surf_spec, n_surf_spec, is_suprathermal_on)
    ALLOCATE(s(total_species), STAT=alloc_stat)
    IF (alloc_stat /= 0) CALL handle_alloc_error('species array', total_species)
    IF (MODEL_EXPERIMENT == 0) THEN
      ALLOCATE(abundances_bulk(total_species, MAX_MONOLAYERS_DEFAULT), STAT=alloc_stat)
      ALLOCATE(timesteps_nml(3000), STAT=alloc_stat)
    ELSE
      ALLOCATE(abundances_bulk(total_species, max_monolayers), STAT=alloc_stat)
      ALLOCATE(timesteps_nml(max_monolayers), STAT=alloc_stat)
    END IF
    IF (alloc_stat /= 0) CALL handle_alloc_error('abundances_bulk or timesteps_nml', total_species)
    DO i = 1, total_species
      ALLOCATE(s(i)%abundance_out(timesteps), STAT=alloc_stat)
      IF (alloc_stat /= 0) CALL handle_alloc_error('abundance_out', i)
      s(i)%abundance_out(:) = 0.0_wp
      s(i)%edes = 0.0_wp
      s(i)%racc = 0.0_wp
      s(i)%rdes = 0.0_wp
      s(i)%abundance = 0.0_wp
      s(i)%frac_abundance = 0.0_wp
    END DO
  END SUBROUTINE allocate_arrays

  SUBROUTINE allocate_reaction_arrays(is_suprathermal_on)
    LOGICAL, INTENT(IN) :: is_suprathermal_on
    INTEGER :: total_reactions, alloc_stat
    total_reactions = nreactions + MERGE(n_surf_react, 0, bulk_chemistry > 0)
    ALLOCATE(r(total_reactions), STAT=alloc_stat)
    IF (alloc_stat /= 0) CALL handle_alloc_error('reaction array', total_reactions)
    ALLOCATE(mre_terms(total_reactions), STAT=alloc_stat)
    ALLOCATE(rd_v2_terms(5, total_reactions), STAT=alloc_stat)
    IF (alloc_stat /= 0) CALL handle_alloc_error('mre_terms or rd_v2_terms', total_reactions)
    rd_v2_terms = 0.0_wp
  END SUBROUTINE allocate_reaction_arrays

  SUBROUTINE link_gas_counterparts()
    INTEGER :: i, j
    DO i = 1, nspecies
      IF (s(i)%name(1:1) == 'g') THEN
        DO j = 1, nspecies
          IF (s(j)%name == s(i)%name(2:LEN_TRIM(s(i)%name))) THEN
            s(i)%gas_idx = j
            s(i)%enthalpia = 0.0_wp
            s(i)%enthalpia_known = 0
          END IF
        END DO
      END IF
    END DO
  END SUBROUTINE link_gas_counterparts

  SUBROUTINE add_suprathermal_surface_species()
    INTEGER :: i
    DO i = 1, n_surf_spec
      s(nspecies + i)%name = TRIM(s(first_surf_spec + i - 1)%name) // '*'
      s(nspecies + i)%idx = nspecies + i
      s(nspecies + i)%gas_idx = s(first_surf_spec + i - 1)%gas_idx
      s(nspecies + i)%weight = s(first_surf_spec + i - 1)%weight
      s(nspecies + i)%natoms = s(first_surf_spec + i - 1)%natoms
      s(nspecies + i)%enthalpia_known = s(first_surf_spec + i - 1)%enthalpia_known
      s(nspecies + i)%enthalpia = s(first_surf_spec + i - 1)%enthalpia
      s(nspecies + i)%edes = s(first_surf_spec + i - 1)%edes
      s(nspecies + i)%racc = s(first_surf_spec + i - 1)%racc
    END DO
    nspecies = nspecies + n_surf_spec
  END SUBROUTINE add_suprathermal_surface_species

  SUBROUTINE add_bulk_species(first_suprathermal_species)
    INTEGER, INTENT(OUT) :: first_suprathermal_species
    INTEGER :: i
    DO i = 1, n_surf_spec
      s(nspecies + i)%name = 'b' // s(first_surf_spec + i - 1)%name(2:LEN_TRIM(s(first_surf_spec + i - 1)%name))
      s(nspecies + i)%idx = nspecies + i
      s(nspecies + i)%gas_idx = s(first_surf_spec + i - 1)%gas_idx
      s(nspecies + i)%weight = s(first_surf_spec + i - 1)%weight
      s(nspecies + i)%natoms = s(first_surf_spec + i - 1)%natoms
      s(nspecies + i)%enthalpia_known = s(first_surf_spec + i - 1)%enthalpia_known
      s(nspecies + i)%enthalpia = s(first_surf_spec + i - 1)%enthalpia
      s(nspecies + i)%edes = s(first_surf_spec + i - 1)%edes
      s(nspecies + i)%racc = s(first_surf_spec + i - 1)%racc
    END DO
    nspecies = nspecies + n_surf_spec
    first_suprathermal_species = nspecies + 1
  END SUBROUTINE add_bulk_species

  SUBROUTINE add_bulk_reactions(output_unit)
    INTEGER, INTENT(INOUT) :: output_unit
    INTEGER :: i, s_r_counter
    IF (bulk_chemistry <= 0) RETURN
    s_r_counter = 0
    DO i = 1, nreactions
      IF (r(i)%r1(1:1) == 'g' .AND. r(i)%r2 /= 'FREEZE' .AND. r(i)%r2 /= 'DESORB' .AND. r(i)%p1(1:1) == 'g') THEN
        s_r_counter = s_r_counter + 1
        CALL copy_reaction_with_bulk_prefix(i, nreactions + s_r_counter)
      END IF
    END DO
    nreactions = nreactions + s_r_counter
    CALL write_reactions_to_file(output_unit, 'bulk_reactions.out', 1, nreactions)
  END SUBROUTINE add_bulk_reactions

  
  SUBROUTINE resize_reaction_array(extra_size)
    INTEGER, INTENT(IN) :: extra_size
    TYPE(reaction), DIMENSION(:), ALLOCATABLE :: rtemp
    INTEGER :: alloc_stat
    ALLOCATE(rtemp(SIZE(r) + extra_size), STAT=alloc_stat)
    IF (alloc_stat /= 0) CALL handle_alloc_error('rtemp for resize', SIZE(r) + extra_size)
    rtemp(1:SIZE(r)) = r
    DEALLOCATE(r)
    ALLOCATE(r(SIZE(rtemp)), STAT=alloc_stat)
    IF (alloc_stat /= 0) CALL handle_alloc_error('r after resize', SIZE(rtemp))
    r = rtemp
    DEALLOCATE(rtemp)
  END SUBROUTINE resize_reaction_array

  SUBROUTINE add_suprathermal_reaction(i, s_r_counter)
    INTEGER, INTENT(IN) :: i
    INTEGER, INTENT(INOUT) :: s_r_counter
    INTEGER :: ii
    IF (r(i)%r1 == r(i)%r2) THEN
      s_r_counter = s_r_counter + 1
      CALL copy_reaction_with_suprathermal(i, nreactions + s_r_counter, 1)
    ELSE
      DO ii = 1, 2
        s_r_counter = s_r_counter + 1
        CALL copy_reaction_with_suprathermal(i, nreactions + s_r_counter, ii)
      END DO
    END IF
  END SUBROUTINE add_suprathermal_reaction

  SUBROUTINE add_reaction_file(input_file, output_file, Nlines, output_unit)
    CHARACTER(LEN=*), INTENT(IN) :: input_file, output_file
    INTEGER, INTENT(OUT) :: Nlines
    INTEGER, INTENT(INOUT) :: output_unit
    INTEGER :: io_stat, file_unit, ii
    LOGICAL :: file_exists
    TYPE(reaction), DIMENSION(:), ALLOCATABLE :: rtemp
    INQUIRE(FILE=input_file, EXIST=file_exists)
    IF (.NOT. file_exists) RETURN
    CALL find_free_unit(file_unit, MIN_UNIT, MAX_UNIT)
    OPEN(file_unit, FILE=input_file, STATUS='OLD', IOSTAT=io_stat)
    IF (io_stat /= 0) RETURN
    Nlines = 0
    DO
      READ(file_unit, *, IOSTAT=io_stat)
      IF (io_stat /= 0) EXIT
      Nlines = Nlines + 1
    END DO
    CLOSE(file_unit)
    IF (Nlines == 0) RETURN
    ALLOCATE(rtemp(nreactions + Nlines))
    rtemp(1:SIZE(r)) = r
    OPEN(file_unit, FILE=input_file, STATUS='OLD', IOSTAT=io_stat)
    DO ii = nreactions + 1, nreactions + Nlines
      READ(file_unit, 1000, IOSTAT=io_stat) rtemp(ii)%idx, rtemp(ii)%r1, rtemp(ii)%r2, &
        rtemp(ii)%p1, rtemp(ii)%p2, rtemp(ii)%p3, rtemp(ii)%p4, rtemp(ii)%p5, &
        rtemp(ii)%alpha, rtemp(ii)%beta, rtemp(ii)%gamma
      IF (io_stat /= 0) CALL handle_io_error('reading ' // TRIM(input_file), ii - nreactions)
      rtemp(ii)%idx = ii
      CALL assign_species_indices_temp(rtemp(ii))
      rtemp(ii)%rtype = get_rtype(rtemp(ii)%r1, rtemp(ii)%r2)
      rtemp(ii)%exothermicity_known = 0
      rtemp(ii)%exothermicity = 0.0_wp
    END DO
    CLOSE(file_unit)
    CALL resize_and_copy_reactions(rtemp)
    CALL write_reactions_to_file(output_unit, output_file, nreactions + 1, nreactions + Nlines)
    nreactions = nreactions + Nlines
1000 FORMAT(1X,I4,1X,2(A10),10X,5(A10),E12.4,1X,F8.2,1X,F8.1)
  END SUBROUTINE add_reaction_file



  SUBROUTINE resize_and_copy_reactions(rtemp)
    TYPE(reaction), DIMENSION(:), INTENT(INOUT), ALLOCATABLE :: rtemp
    INTEGER :: alloc_stat
    DEALLOCATE(r)
    ALLOCATE(r(SIZE(rtemp)), STAT=alloc_stat)
    IF (alloc_stat /= 0) CALL handle_alloc_error('r after resize', SIZE(rtemp))
    r = rtemp
    DEALLOCATE(rtemp)
  END SUBROUTINE resize_and_copy_reactions


  SUBROUTINE assign_species_indices_temp(reac)
    TYPE(reaction), INTENT(INOUT) :: reac
    reac%ir1 = species_idx_with_autocreate(reac%r1)  ! Use the new function
    reac%ir2 = species_idx_with_autocreate(reac%r2)  ! Use the new function
    reac%ip1 = species_idx_with_autocreate(reac%p1)  ! Use the new function
    reac%ip2 = species_idx_with_autocreate(reac%p2)  ! Use the new function
    reac%ip3 = species_idx_with_autocreate(reac%p3)  ! Use the new function
    reac%ip4 = species_idx_with_autocreate(reac%p4)  ! Use the new function
    reac%ip5 = species_idx_with_autocreate(reac%p5)  ! Use the new function
    IF (ANY([reac%ir1, reac%ir2, reac%ip1, reac%ip2, reac%ip3, reac%ip4, reac%ip5] == -1)) THEN
      WRITE(*,'(A,I5,1X,7A10)') 'ERROR: Unknown species in reaction: ', reac%idx, &
        reac%r1, reac%r2, reac%p1, reac%p2, reac%p3, reac%p4, reac%p5
      CALL cleanup_and_exit()
    END IF
  END SUBROUTINE assign_species_indices_temp

  SUBROUTINE copy_reaction_with_bulk_prefix(src_idx, dest_idx)
    INTEGER, INTENT(IN) :: src_idx, dest_idx
    r(dest_idx) = r(src_idx)
    r(dest_idx)%idx = dest_idx
    IF (r(src_idx)%r1(1:1) == 'g') r(dest_idx)%r1 = 'b' // r(src_idx)%r1(2:LEN_TRIM(r(src_idx)%r1))
    IF (r(src_idx)%r2(1:1) == 'g') r(dest_idx)%r2 = 'b' // r(src_idx)%r2(2:LEN_TRIM(r(src_idx)%r2))
    IF (r(src_idx)%p1(1:1) == 'g') r(dest_idx)%p1 = 'b' // r(src_idx)%p1(2:LEN_TRIM(r(src_idx)%p1))
    IF (r(src_idx)%p2(1:1) == 'g') r(dest_idx)%p2 = 'b' // r(src_idx)%p2(2:LEN_TRIM(r(src_idx)%p2))
    IF (r(src_idx)%p3(1:1) == 'g') r(dest_idx)%p3 = 'b' // r(src_idx)%p3(2:LEN_TRIM(r(src_idx)%p3))
    IF (r(src_idx)%p4(1:1) == 'g') r(dest_idx)%p4 = 'b' // r(src_idx)%p4(2:LEN_TRIM(r(src_idx)%p4))
    IF (r(src_idx)%p5(1:1) == 'g') r(dest_idx)%p5 = 'b' // r(src_idx)%p5(2:LEN_TRIM(r(src_idx)%p5))
    CALL assign_species_indices(dest_idx)
    r(dest_idx)%rtype = get_rtype(r(dest_idx)%r1, r(dest_idx)%r2)
  END SUBROUTINE copy_reaction_with_bulk_prefix

  SUBROUTINE copy_reaction_with_suprathermal(src_idx, dest_idx, reactant_num)
    INTEGER, INTENT(IN) :: src_idx, dest_idx, reactant_num
    r(dest_idx) = r(src_idx)
    r(dest_idx)%idx = dest_idx
    IF (reactant_num == 1) r(dest_idx)%r1 = TRIM(r(src_idx)%r1) // '*'
    IF (reactant_num == 2) r(dest_idx)%r2 = TRIM(r(src_idx)%r2) // '*'
    CALL assign_species_indices(dest_idx)
    r(dest_idx)%rtype = get_rtype(r(dest_idx)%r1, r(dest_idx)%r2)
    r(dest_idx)%exothermicity_known = 0
    r(dest_idx)%exothermicity = 0.0_wp
  END SUBROUTINE copy_reaction_with_suprathermal


  SUBROUTINE write_species_output(file_unit)
    INTEGER, INTENT(INOUT) :: file_unit
    INTEGER :: i, io_stat
    CALL find_free_unit(file_unit, MIN_UNIT, MAX_UNIT)
    OPEN(file_unit, FILE='species.out', STATUS='REPLACE', IOSTAT=io_stat)
    IF (io_stat == 0) THEN
      WRITE(file_unit, '(A)') 'Species, Weight, Number of Atoms, Index, Gas twin index, Desorption energy'
      DO i = 1, nspecies
        WRITE(file_unit, '(A10,1X,1pE12.4,1X,I3,1X,2I4,1pE12.4)') &
          s(i)%name, s(i)%weight, s(i)%natoms, s(i)%idx, s(i)%gas_idx, s(i)%edes
      END DO
      CLOSE(file_unit)
    END IF
  END SUBROUTINE write_species_output

  SUBROUTINE verify_reactions()
    INTEGER :: i, prodatoms, reactatoms
    CHARACTER(LEN=12), PARAMETER :: special_reactions(*) = &
      ['QUENCH', 'CRPHOT', 'PHOTON', 'FREEZE', 'DESORB', 'IONRAD', &
       'G-    ', 'G0    ', 'CR    ', 'CRP   ', 'PHOION', 'PHOEXC']
    DO i = 1, nreactions
      IF (r(i)%ir1 == 0) CYCLE
      prodatoms = s(r(i)%ip1)%natoms
      IF (r(i)%ip2 /= 0) prodatoms = prodatoms + s(r(i)%ip2)%natoms
      IF (r(i)%ip3 /= 0) prodatoms = prodatoms + s(r(i)%ip3)%natoms
      IF (r(i)%ip4 /= 0) prodatoms = prodatoms + s(r(i)%ip4)%natoms
      IF (r(i)%ip5 /= 0) prodatoms = prodatoms + s(r(i)%ip5)%natoms
      reactatoms = s(r(i)%ir1)%natoms
      IF (r(i)%ir2 /= 0 .AND. .NOT. ANY(special_reactions == r(i)%r2)) &
        reactatoms = reactatoms + s(r(i)%ir2)%natoms
      IF (prodatoms /= reactatoms .AND. r(i)%r1(1:2) /= 'e-') THEN
        WRITE(*,'(A,I0,A,I0)') 'Debug: Reaction ', r(i)%idx, ' reactatoms=', reactatoms
        WRITE(*,'(A,I0)') 'Debug: prodatoms=', prodatoms
        WRITE(*,'(A,I0,A,I0,A,I0)') 'WARNING: Atom imbalance in reaction ', r(i)%idx, &
          ': Reactants=', reactatoms, ' Products=', prodatoms
      END IF
      IF (ISNAN(r(i)%alpha)) r(i)%alpha = 0.0_wp
      IF (ISNAN(r(i)%beta)) r(i)%beta = 0.0_wp
      IF (ISNAN(r(i)%gamma)) r(i)%gamma = 0.0_wp
      IF (ISNAN(r(i)%rate)) r(i)%rate = 0.0_wp
      IF (ISNAN(r(i)%exothermicity)) r(i)%exothermicity = 0.0_wp
    END DO
  END SUBROUTINE verify_reactions

  SUBROUTINE handle_io_error(context, position)
    CHARACTER(LEN=*), INTENT(IN) :: context
    INTEGER, INTENT(IN), OPTIONAL :: position
    IF (PRESENT(position)) THEN
      WRITE(*,'(A,A,I0)') 'ERROR: Failed ', TRIM(context), ' at position ', position
    ELSE
      WRITE(*,'(A,A)') 'ERROR: Failed ', TRIM(context)
    END IF
    CALL cleanup_and_exit()
  END SUBROUTINE handle_io_error

  SUBROUTINE handle_alloc_error(array_name, size_val)
    CHARACTER(LEN=*), INTENT(IN) :: array_name
    INTEGER, INTENT(IN) :: size_val
    WRITE(*,'(A,A,A,I0)') 'ERROR: Failed to allocate ', TRIM(array_name), ' with size ', size_val
    CALL cleanup_and_exit()
  END SUBROUTINE handle_alloc_error

  SUBROUTINE cleanup_and_exit()
    INTEGER :: i
    IF (ALLOCATED(s)) THEN
      DO i = 1, SIZE(s)
        IF (ALLOCATED(s(i)%abundance_out)) DEALLOCATE(s(i)%abundance_out)
      END DO
      DEALLOCATE(s)
    END IF
    IF (ALLOCATED(abundances_bulk)) DEALLOCATE(abundances_bulk)
    IF (ALLOCATED(timesteps_nml)) DEALLOCATE(timesteps_nml)
    IF (ALLOCATED(r)) DEALLOCATE(r)
    IF (ALLOCATED(mre_terms)) DEALLOCATE(mre_terms)
    IF (ALLOCATED(rd_v2_terms)) DEALLOCATE(rd_v2_terms)
    STOP
  END SUBROUTINE cleanup_and_exit

  SUBROUTINE read_enthalpias
    IMPLICIT NONE
    CHARACTER(LEN=10) :: local_species
    REAL(KIND=wp) :: local_enthalpia
    INTEGER :: local_enthalpia_known, i, io_stat, unit_num
    LOGICAL :: file_exists
    INQUIRE(FILE='enthalpias.txt', EXIST=file_exists)
    IF (.NOT. file_exists) THEN
      WRITE(*,*) 'WARNING: enthalpias.txt not found'
      RETURN
    END IF
    CALL find_free_unit(unit_num, MIN_UNIT, MAX_UNIT)
    OPEN(unit_num, FILE='enthalpias.txt', STATUS='OLD', IOSTAT=io_stat)
    IF (io_stat /= 0) RETURN
    READ(unit_num, *, IOSTAT=io_stat) ! Skip header
    IF (io_stat /= 0) THEN
      CLOSE(unit_num)
      RETURN
    END IF
    DO
      READ(unit_num, '(A10,1X,E12.4,2X,I1)', IOSTAT=io_stat) &
        local_species, local_enthalpia, local_enthalpia_known
      IF (io_stat /= 0) EXIT
      DO i = 1, nspecies
        IF (s(i)%name == 'g' // TRIM(local_species) .OR. s(i)%name == TRIM(local_species)) THEN
          s(i)%enthalpia = local_enthalpia * 1.0e3_wp / 8.31_wp
          s(i)%enthalpia_known = local_enthalpia_known
        END IF
      END DO
    END DO
    CLOSE(unit_num)
  END SUBROUTINE read_enthalpias

  SUBROUTINE get_reaction_thermodynamics(ii)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ii
    IF (ii <= 0 .OR. ii > SIZE(r)) RETURN
    r(ii)%exothermicity_known = s(r(ii)%ir1)%enthalpia_known
    IF (r(ii)%ir2 /= 0) r(ii)%exothermicity_known = r(ii)%exothermicity_known * s(r(ii)%ir2)%enthalpia_known
    IF (r(ii)%ip1 /= 0) r(ii)%exothermicity_known = r(ii)%exothermicity_known * s(r(ii)%ip1)%enthalpia_known
    IF (r(ii)%ip2 /= 0) r(ii)%exothermicity_known = r(ii)%exothermicity_known * s(r(ii)%ip2)%enthalpia_known
    IF (r(ii)%ip3 /= 0) r(ii)%exothermicity_known = r(ii)%exothermicity_known * s(r(ii)%ip3)%enthalpia_known
    IF (r(ii)%ip4 /= 0) r(ii)%exothermicity_known = r(ii)%exothermicity_known * s(r(ii)%ip4)%enthalpia_known
    IF (r(ii)%ip5 /= 0) r(ii)%exothermicity_known = r(ii)%exothermicity_known * s(r(ii)%ip5)%enthalpia_known
    r(ii)%exothermicity = 0.0_wp
    IF (r(ii)%exothermicity_known == 1) THEN
      r(ii)%exothermicity = r(ii)%exothermicity - s(r(ii)%ir1)%enthalpia
      IF (r(ii)%ir2 /= 0) r(ii)%exothermicity = r(ii)%exothermicity - s(r(ii)%ir2)%enthalpia
      IF (r(ii)%ip1 /= 0) r(ii)%exothermicity = r(ii)%exothermicity + s(r(ii)%ip1)%enthalpia
      IF (r(ii)%ip2 /= 0) r(ii)%exothermicity = r(ii)%exothermicity + s(r(ii)%ip2)%enthalpia
      IF (r(ii)%ip3 /= 0) r(ii)%exothermicity = r(ii)%exothermicity + s(r(ii)%ip3)%enthalpia
      IF (r(ii)%ip4 /= 0) r(ii)%exothermicity = r(ii)%exothermicity + s(r(ii)%ip4)%enthalpia
      IF (r(ii)%ip5 /= 0) r(ii)%exothermicity = r(ii)%exothermicity + s(r(ii)%ip5)%enthalpia
      r(ii)%exothermicity = -r(ii)%exothermicity
    END IF
  END SUBROUTINE get_reaction_thermodynamics

  SUBROUTINE get_rd_efficiency(ii)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ii
    REAL(KIND=wp) :: P
    IF (ii <= 0 .OR. ii > SIZE(r)) RETURN
    SELECT CASE (des_reactive_type)
    CASE (0) ! Garrod_ea07
      IF (r(ii)%rtype == 13) THEN
        IF (r(ii)%p1(1:1) == 'g') THEN
          IF (r(ii)%ip2 /= 0) THEN
            r(ii)%alpha = 1.0_wp * r(ii)%alpha
          ELSE
            IF (r(ii)%exothermicity_known == 1 .AND. r(ii)%ip1 > 0) THEN
              P = (1.0_wp - s(r(ii)%ip1)%edes / r(ii)%exothermicity) ** &
                  (MAX(2.0_wp, 3.0_wp * s(r(ii)%ip1)%natoms - 5.0_wp) - 1.0_wp)
              IF (P < 0.0_wp .OR. ISNAN(P)) P = 0.0_wp
              IF (P > 1.0_wp) P = 1.0_wp
              r(ii)%alpha = (1.0_wp - (des_reactive * P) / (1.0_wp + des_reactive * P)) * r(ii)%alpha
            ELSE
              r(ii)%alpha = (1.0_wp - des_reactive) * r(ii)%alpha
            END IF
          END IF
        ELSE
          IF (r(ii)%ip2 /= 0) THEN
            r(ii)%alpha = 0.0_wp
          ELSE
            IF (r(ii)%exothermicity_known == 1 .AND. r(ii)%ip1 > 0) THEN
              P = (1.0_wp - s(r(ii)%ip1)%edes / r(ii)%exothermicity) ** &
                  (MAX(2.0_wp, 3.0_wp * s(r(ii)%ip1)%natoms - 5.0_wp) - 1.0_wp)
              IF (P < 0.0_wp .OR. ISNAN(P)) P = 0.0_wp
              IF (P > 1.0_wp) P = 1.0_wp
              r(ii)%alpha = (des_reactive * P) / (1.0_wp + des_reactive * P) * r(ii)%alpha
            ELSE
              r(ii)%alpha = des_reactive * r(ii)%alpha
            END IF
          END IF
        END IF
      END IF
    CASE (1, 2) ! Vasyunin&Herbst13 or Minissale&Dulieu
      IF (r(ii)%rtype == 13) THEN
        IF (r(ii)%p1(1:1) == 'g') THEN
          r(ii)%alpha = (1.0_wp - des_reactive) * r(ii)%alpha
        ELSE
          r(ii)%alpha = des_reactive * r(ii)%alpha
        END IF
      END IF
    CASE DEFAULT
      WRITE(*,'(A,I0)') 'ERROR: Unknown des_reactive_type: ', des_reactive_type
      CALL cleanup_and_exit()
    END SELECT
  END SUBROUTINE get_rd_efficiency

  INTEGER FUNCTION get_rtype(r1, r2)
    IMPLICIT NONE
    CHARACTER(LEN=10), INTENT(IN) :: r1, r2
    LOGICAL :: r1IsIon, r2IsIon
    INTEGER :: r1_len, r2_len
    get_rtype = 1
    r1_len = LEN_TRIM(r1)
    r2_len = LEN_TRIM(r2)
    IF (r1_len == 0 .OR. r2_len == 0) RETURN
    IF (r2 == 'CRP') get_rtype = 2
    IF (r2 == 'PHOTON') THEN
      IF (r1 == 'H2') THEN
        get_rtype = 4
      ELSE IF (r1 == 'CO') THEN
        get_rtype = 5
      ELSE
        get_rtype = 3
      END IF
    END IF
    IF (r2 == 'CRPHOT') get_rtype = 6
    IF (r2 == 'G-') get_rtype = 7
    IF (r2 == 'G0') get_rtype = 8
    IF (r1 == 'G0') get_rtype = 9
    IF (r1 == 'G+') get_rtype = 10
    IF (r2 == 'FREEZE') get_rtype = 11
    IF (r2 == 'DESORB') get_rtype = 12
    IF (r1(1:1) == 'g' .AND. r2(1:1) == 'g') get_rtype = 13
    IF (r1(1:1) == 'b' .AND. r2(1:1) == 'b') get_rtype = 14
    IF ((r1(r1_len:r1_len) == '*' .OR. r2(r2_len:r2_len) == '*') .AND. &
        (r1(1:1) == 'g' .OR. r2(1:1) == 'g')) get_rtype = 15
    IF ((r1(r1_len:r1_len) == '*' .OR. r2(r2_len:r2_len) == '*') .AND. &
        (r1(1:1) == 'b' .OR. r2(1:1) == 'b')) get_rtype = 16
    IF (r2 == 'IONRAD') get_rtype = 17
    IF (r2 == 'QUENCH') get_rtype = 18
    IF (r2 == 'PHOION') get_rtype = 19
    IF (r2 == 'PHOEXC') get_rtype = 20
    r1IsIon = (r1(r1_len:r1_len) == '+' .OR. r1(r1_len:r1_len) == '-')
    r2IsIon = (r2(r2_len:r2_len) == '+' .OR. r2(r2_len:r2_len) == '-')
    IF ((r1IsIon .AND. .NOT. r2IsIon) .OR. (.NOT. r1IsIon .AND. r2IsIon)) THEN
      IF (r1(1:1) == 'g') get_rtype = 21
      IF (r1(1:1) == 'b') get_rtype = 22
    END IF
    IF (r1IsIon .AND. r2IsIon) THEN
      IF (r1(1:1) == 'g') get_rtype = 23
      IF (r1(1:1) == 'b') get_rtype = 24
    END IF
  END FUNCTION get_rtype

  REAL(KIND=wp) FUNCTION REALVALUE(str)
    CHARACTER(LEN=*), INTENT(IN) :: str
    INTEGER :: io_stat
    READ(str, *, IOSTAT=io_stat) REALVALUE
    IF (io_stat /= 0) REALVALUE = 0.0_wp
  END FUNCTION REALVALUE



  SUBROUTINE extend_species_array(extra_slots)
    INTEGER, INTENT(IN) :: extra_slots
    TYPE(species), DIMENSION(:), ALLOCATABLE :: stemp
    INTEGER :: i, j, k, alloc_stat, old_size, new_size
    
    old_size = SIZE(s)
    new_size = old_size + extra_slots
    
    ALLOCATE(stemp(new_size), STAT=alloc_stat)
    IF (alloc_stat /= 0) CALL handle_alloc_error('stemp for extend', new_size)
    
    ! Copy existing data (except abundance_out)
    DO i = 1, old_size
      ! Copy all members except abundance_out
      stemp(i)%name = s(i)%name
      stemp(i)%idx = s(i)%idx
      stemp(i)%gas_idx = s(i)%gas_idx
      stemp(i)%weight = s(i)%weight
      stemp(i)%natoms = s(i)%natoms
      stemp(i)%enthalpia_known = s(i)%enthalpia_known
      stemp(i)%enthalpia = s(i)%enthalpia
      stemp(i)%edes = s(i)%edes
      stemp(i)%racc = s(i)%racc
      stemp(i)%rdes = s(i)%rdes
      stemp(i)%abundance = s(i)%abundance
      stemp(i)%frac_abundance = s(i)%frac_abundance
      
      ! Now allocate and copy abundance_out
      ALLOCATE(stemp(i)%abundance_out(timesteps), STAT=alloc_stat)
      IF (alloc_stat /= 0) CALL handle_alloc_error('abundance_out in extend', i)
      stemp(i)%abundance_out(:) = s(i)%abundance_out(:)
    END DO
    
    ! Initialize new entries - Use a different loop variable (j)
    DO j = old_size + 1, new_size
      stemp(j)%name = ' '
      stemp(j)%idx = 0
      stemp(j)%gas_idx = 0
      stemp(j)%weight = 0.0_wp
      stemp(j)%natoms = 0
      stemp(j)%enthalpia_known = 0
      stemp(j)%enthalpia = 0.0_wp
      stemp(j)%edes = 0.0_wp
      stemp(j)%racc = 0.0_wp
      stemp(j)%rdes = 0.0_wp
      stemp(j)%abundance = 0.0_wp
      stemp(j)%frac_abundance = 0.0_wp
      ALLOCATE(stemp(j)%abundance_out(timesteps), STAT=alloc_stat)
      IF (alloc_stat /= 0) CALL handle_alloc_error('abundance_out in extend', j)
      stemp(j)%abundance_out(:) = 0.0_wp
    END DO
    
    ! Replace the old array with the new one
    DO i = 1, old_size
      DEALLOCATE(s(i)%abundance_out)
    END DO
    DEALLOCATE(s)
    ALLOCATE(s(new_size), STAT=alloc_stat)
    IF (alloc_stat /= 0) CALL handle_alloc_error('s after extend', new_size)
    
    ! Copy from temp array to new s array - Use a different loop variable (k)
    DO k = 1, new_size
      s(k)%name = stemp(k)%name
      s(k)%idx = stemp(k)%idx
      s(k)%gas_idx = stemp(k)%gas_idx
      s(k)%weight = stemp(k)%weight
      s(k)%natoms = stemp(k)%natoms
      s(k)%enthalpia_known = stemp(k)%enthalpia_known
      s(k)%enthalpia = stemp(k)%enthalpia
      s(k)%edes = stemp(k)%edes
      s(k)%racc = stemp(k)%racc
      s(k)%rdes = stemp(k)%rdes
      s(k)%abundance = stemp(k)%abundance
      s(k)%frac_abundance = stemp(k)%frac_abundance
      
      ALLOCATE(s(k)%abundance_out(timesteps), STAT=alloc_stat)
      IF (alloc_stat /= 0) CALL handle_alloc_error('abundance_out in new s', k)
      s(k)%abundance_out(:) = stemp(k)%abundance_out(:)
      
      DEALLOCATE(stemp(k)%abundance_out)
    END DO
    DEALLOCATE(stemp)
    
    WRITE(*,'(A,I0,A,I0)') 'Extended species array from ', old_size, ' to ', new_size
  END SUBROUTINE extend_species_array

  SUBROUTINE add_quenching_reactions(file_unit, output_unit)
    INTEGER, INTENT(INOUT) :: file_unit, output_unit
    INTEGER :: ii, jj, io_stat, first_suprathermal_species
    CHARACTER(LEN=10) :: groundstate
    TYPE(reaction), DIMENSION(:), ALLOCATABLE :: rtemp
    INTEGER :: quench_count, actual_unit
    
    ! Find the starting index for suprathermal species
    first_suprathermal_species = 0
    DO ii = 1, nspecies
      IF (LEN_TRIM(s(ii)%name) > 1) THEN
        IF (s(ii)%name(LEN_TRIM(s(ii)%name):LEN_TRIM(s(ii)%name)) == '*') THEN
          first_suprathermal_species = ii
          EXIT
        END IF
      END IF
    END DO
    
    ! Check if we found any suprathermal species
    IF (first_suprathermal_species <= 0) THEN
      WRITE(*,*) 'WARNING: No suprathermal species found. Skipping quenching reactions.'
      RETURN
    END IF
    
    ! Calculate number of quenching reactions
    quench_count = 0
    DO ii = first_suprathermal_species, nspecies
      IF (LEN_TRIM(s(ii)%name) > 1) THEN
        IF (s(ii)%name(LEN_TRIM(s(ii)%name):LEN_TRIM(s(ii)%name)) == '*') THEN
          quench_count = quench_count + 1
        END IF
      END IF
    END DO
    
    IF (quench_count <= 0) RETURN
    
    ALLOCATE(rtemp(nreactions + quench_count))
    rtemp(1:SIZE(r)) = r
    
    ! Ensure we have a valid unit number
    actual_unit = MIN_UNIT  ! Default to a valid unit
    CALL find_free_unit(actual_unit, MIN_UNIT, MAX_UNIT)
    IF (actual_unit == -1) THEN
      WRITE(*,*) 'WARNING: Could not find free unit for quenching.out. Using unit ', MIN_UNIT+1
      actual_unit = MIN_UNIT + 1
    END IF
    
    ! Open output file
    OPEN(UNIT=actual_unit, FILE='quenching.out', STATUS='REPLACE', IOSTAT=io_stat)
    IF (io_stat /= 0) THEN
      WRITE(*,'(A,I0)') 'WARNING: Failed to open quenching.out file, IOSTAT=', io_stat
      RETURN
    END IF
    
    ! Create quenching reactions
    jj = 0
    DO ii = first_suprathermal_species, nspecies
      IF (LEN_TRIM(s(ii)%name) > 1) THEN
        IF (s(ii)%name(LEN_TRIM(s(ii)%name):LEN_TRIM(s(ii)%name)) == '*') THEN
          jj = jj + 1
          IF (jj > quench_count) EXIT  ! Safety check
          
          groundstate = s(ii)%name(1:LEN_TRIM(s(ii)%name) - 1)
          rtemp(nreactions + jj)%r1 = s(ii)%name
          rtemp(nreactions + jj)%r2 = 'QUENCH'
          rtemp(nreactions + jj)%p1 = groundstate
          rtemp(nreactions + jj)%p2 = ' '
          rtemp(nreactions + jj)%p3 = ' '
          rtemp(nreactions + jj)%p4 = ' '
          rtemp(nreactions + jj)%p5 = ' '
          rtemp(nreactions + jj)%idx = nreactions + jj
          
          ! Use normal lookup first to avoid auto-creation for reactants
          rtemp(nreactions + jj)%ir1 = species_idx(rtemp(nreactions + jj)%r1)
          rtemp(nreactions + jj)%ir2 = 0
          rtemp(nreactions + jj)%ip1 = species_idx(groundstate)
          rtemp(nreactions + jj)%ip2 = 0
          rtemp(nreactions + jj)%ip3 = 0
          rtemp(nreactions + jj)%ip4 = 0
          rtemp(nreactions + jj)%ip5 = 0
          
          ! If we couldn't find the species, try with auto-creation
          IF (rtemp(nreactions + jj)%ir1 == -1) THEN
            rtemp(nreactions + jj)%ir1 = species_idx_with_autocreate(rtemp(nreactions + jj)%r1)
          END IF
          IF (rtemp(nreactions + jj)%ip1 == -1) THEN
            rtemp(nreactions + jj)%ip1 = species_idx_with_autocreate(groundstate)
          END IF
          
          ! Skip if we still can't find the species
          IF (rtemp(nreactions + jj)%ir1 == -1 .OR. rtemp(nreactions + jj)%ip1 == -1) THEN
            WRITE(*,'(A,I0,A,A10,A,A10)') 'WARNING: Skipping quenching reaction for ', &
              nreactions + jj, ', r1=', rtemp(nreactions + jj)%r1, ', p1=', groundstate
            CYCLE
          END IF
          
          rtemp(nreactions + jj)%alpha = 1.0_wp
          rtemp(nreactions + jj)%beta = 1.0_wp
          rtemp(nreactions + jj)%gamma = 1.0_wp
          rtemp(nreactions + jj)%rtype = get_rtype(rtemp(nreactions + jj)%r1, rtemp(nreactions + jj)%r2)
          rtemp(nreactions + jj)%exothermicity = 0.0_wp
          rtemp(nreactions + jj)%exothermicity_known = 0
          
          ! FIX: Use the exact same FORMAT statement as used elsewhere in the code (FORMAT 1001)
          ! The original code is using FORMAT 1001 in read_rate06database, so we use the same format here
          WRITE(actual_unit, '(1X,I4,1X,I2,1X,2(A10),10X,5(A10),1pE14.5,1X,I1,1X,1pE14.5)', IOSTAT=io_stat) &
            rtemp(nreactions + jj)%idx, rtemp(nreactions + jj)%rtype, &
            rtemp(nreactions + jj)%r1, rtemp(nreactions + jj)%r2, &
            rtemp(nreactions + jj)%p1, rtemp(nreactions + jj)%p2, &
            rtemp(nreactions + jj)%p3, rtemp(nreactions + jj)%p4, &
            rtemp(nreactions + jj)%p5, rtemp(nreactions + jj)%exothermicity, &
            rtemp(nreactions + jj)%exothermicity_known, rtemp(nreactions + jj)%alpha
          
          IF (io_stat /= 0) THEN
            WRITE(*,'(A,I0)') 'WARNING: Failed to write quenching reaction, IOSTAT=', io_stat
          END IF
        END IF
      END IF
    END DO
    
    CLOSE(actual_unit)
    CALL resize_and_copy_reactions(rtemp)
    nreactions = nreactions + jj  ! Only add the successful reactions
    file_unit = actual_unit
  END SUBROUTINE add_quenching_reactions

INTEGER FUNCTION species_idx_with_autocreate(species_name)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: species_name
  INTEGER :: i, j, k, new_idx
  LOGICAL :: is_suprathermal
  CHARACTER(LEN=10) :: base_name
  INTEGER :: alloc_stat, old_size, new_size
  TYPE(species), DIMENSION(:), ALLOCATABLE :: stemp
  
  ! Skip empty species names
  IF (LEN_TRIM(species_name) == 0) THEN
    species_idx_with_autocreate = 0  ! Return 0 for empty species
    RETURN
  END IF
  
  ! First try normal lookup
  species_idx_with_autocreate = species_idx(species_name)
  
  ! If not found and it's a suprathermal bulk species, create it
  IF (species_idx_with_autocreate == -1 .AND. LEN_TRIM(species_name) > 1) THEN
    is_suprathermal = (species_name(1:1) == 'b' .AND. &
                     species_name(LEN_TRIM(species_name):LEN_TRIM(species_name)) == '*')
    
    IF (is_suprathermal) THEN
      ! Extract base species name (without *)
      base_name = species_name(1:LEN_TRIM(species_name)-1)
      
      ! Find the base species
      DO i = 1, nspecies
        IF (s(i)%name == base_name) THEN
          ! Create new species
          new_idx = nspecies + 1
          
          ! Allocate more space if needed
          IF (new_idx > SIZE(s)) THEN
            WRITE(*,'(A)') 'WARNING: Auto-extending species array'
            
            ! FIX: Safer array extension
            old_size = SIZE(s)
            new_size = old_size + 10  ! Add 10 more slots
            
            ! Allocate temporary array
            ALLOCATE(stemp(new_size), STAT=alloc_stat)
            IF (alloc_stat /= 0) THEN
              WRITE(*,'(A)') 'ERROR: Failed to allocate temporary species array'
              species_idx_with_autocreate = -1
              RETURN
            END IF
            
            ! Copy existing data
            stemp(1:old_size) = s(1:old_size)
            
            ! Initialize new entries - Use a different loop variable (j)
            DO j = old_size + 1, new_size
              stemp(j)%name = ' '
              stemp(j)%idx = 0
              stemp(j)%gas_idx = 0
              stemp(j)%weight = 0.0_wp
              stemp(j)%natoms = 0
              stemp(j)%enthalpia_known = 0
              stemp(j)%enthalpia = 0.0_wp
              stemp(j)%edes = 0.0_wp
              stemp(j)%racc = 0.0_wp
              stemp(j)%rdes = 0.0_wp
              stemp(j)%abundance = 0.0_wp
              stemp(j)%frac_abundance = 0.0_wp
              
              ! Allocate abundance_out for each new entry
              ALLOCATE(stemp(j)%abundance_out(timesteps), STAT=alloc_stat)
              IF (alloc_stat /= 0) THEN
                WRITE(*,'(A,I0)') 'ERROR: Failed to allocate abundance_out for new species ', j
                species_idx_with_autocreate = -1
                RETURN
              END IF
              stemp(j)%abundance_out = 0.0_wp
            END DO
            
            ! Replace old array with new one
            DEALLOCATE(s)
            ALLOCATE(s(new_size), STAT=alloc_stat)
            IF (alloc_stat /= 0) THEN
              WRITE(*,'(A)') 'ERROR: Failed to allocate new species array'
              species_idx_with_autocreate = -1
              RETURN
            END IF
            
            s = stemp
            
            ! Clean up temp array (but preserve the abundance_out arrays)
            ! Use a different loop variable (k)
            DO k = 1, new_size
              IF (ALLOCATED(stemp(k)%abundance_out)) THEN
                DEALLOCATE(stemp(k)%abundance_out)
              END IF
            END DO
            DEALLOCATE(stemp)
            
            WRITE(*,'(A,I0,A,I0)') 'Extended species array from ', old_size, ' to ', new_size
          END IF
          
          ! Add the new species
          s(new_idx)%name = species_name
          s(new_idx)%idx = new_idx
          s(new_idx)%gas_idx = s(i)%gas_idx
          s(new_idx)%weight = s(i)%weight
          s(new_idx)%natoms = s(i)%natoms
          s(new_idx)%enthalpia_known = s(i)%enthalpia_known
          s(new_idx)%enthalpia = s(i)%enthalpia
          s(new_idx)%edes = s(i)%edes
          s(new_idx)%racc = s(i)%racc
          s(new_idx)%rdes = 0.0_wp
          s(new_idx)%abundance = 0.0_wp
          s(new_idx)%frac_abundance = 0.0_wp
          
          ! FIX: Check if abundance_out is already allocated
          IF (ALLOCATED(s(new_idx)%abundance_out)) THEN
            DEALLOCATE(s(new_idx)%abundance_out)
          END IF
          
          ! Allocate abundance array with proper error checking
          ALLOCATE(s(new_idx)%abundance_out(timesteps), STAT=alloc_stat)
          IF (alloc_stat /= 0) THEN
            WRITE(*,'(A,I0)') 'ERROR: Failed to allocate abundance_out for species ', new_idx
            species_idx_with_autocreate = -1
            RETURN
          END IF
          s(new_idx)%abundance_out = 0.0_wp
          
          nspecies = new_idx
          species_idx_with_autocreate = new_idx
          
          WRITE(*,'(A,A10,A,I0)') 'AUTO-CREATED missing species: ', species_name, &
                                 ' with index=', new_idx
          EXIT
        END IF
      END DO
    END IF
  END IF
END FUNCTION species_idx_with_autocreate

SUBROUTINE add_suprathermal_bulk_and_reactions(Nsup_g, first_suprathermal_react, output_unit)
  INTEGER, INTENT(OUT) :: Nsup_g, first_suprathermal_react
  INTEGER, INTENT(INOUT) :: output_unit
  INTEGER :: i, ii, j, s_r_counter, Nlines, io_stat, file_unit
  INTEGER :: species_added, base_idx
  CHARACTER(LEN=10) :: base_name, suprathermal_name
  LOGICAL :: species_exists, already_has_asterisk
  TYPE(reaction), DIMENSION(:), ALLOCATABLE :: rtemp
  
  ! First, add suprathermal bulk species for ALL existing bulk species
  species_added = 0
  
  ! Create a mapping from normal bulk species to their suprathermal versions
  WRITE(*,'(A)') 'Creating suprathermal bulk species for all existing bulk species'
  
  ! Start from the first bulk species
  DO i = 1, nspecies
    IF (s(i)%name(1:1) == 'b') THEN
      base_name = s(i)%name
      
      ! Check if the species already has an asterisk
      already_has_asterisk = .FALSE.
      IF (LEN_TRIM(base_name) > 0) THEN
        already_has_asterisk = (base_name(LEN_TRIM(base_name):LEN_TRIM(base_name)) == '*')
      END IF
      
      ! Skip species that already have an asterisk
      IF (.NOT. already_has_asterisk) THEN
        suprathermal_name = TRIM(base_name) // '*'
        
        ! Check if this species already exists
        species_exists = .FALSE.
        DO j = 1, nspecies
          IF (s(j)%name == suprathermal_name) THEN
            species_exists = .TRUE.
            EXIT
          END IF
        END DO
        
        IF (.NOT. species_exists) THEN
          species_added = species_added + 1
          ! FIX: Ensure we have enough space in the species array
          IF (nspecies + species_added > SIZE(s)) THEN
            WRITE(*,'(A)') 'ERROR: Not enough space in species array. Extend it first.'
            CALL extend_species_array(10)  ! Add more slots
          END IF
          
          s(nspecies + species_added)%name = suprathermal_name
          s(nspecies + species_added)%idx = nspecies + species_added
          s(nspecies + species_added)%gas_idx = s(i)%gas_idx
          s(nspecies + species_added)%weight = s(i)%weight
          s(nspecies + species_added)%natoms = s(i)%natoms
          s(nspecies + species_added)%enthalpia_known = s(i)%enthalpia_known
          s(nspecies + species_added)%enthalpia = s(i)%enthalpia
          s(nspecies + species_added)%edes = s(i)%edes
          s(nspecies + species_added)%racc = s(i)%racc
          s(nspecies + species_added)%rdes = 0.0_wp
          s(nspecies + species_added)%abundance = 0.0_wp
          s(nspecies + species_added)%frac_abundance = 0.0_wp
          
          ! FIX: Properly allocate abundance_out
          IF (ALLOCATED(s(nspecies + species_added)%abundance_out)) THEN
            DEALLOCATE(s(nspecies + species_added)%abundance_out)
          END IF
          ALLOCATE(s(nspecies + species_added)%abundance_out(timesteps), STAT=io_stat)
          IF (io_stat /= 0) THEN
            WRITE(*,'(A,I0)') 'ERROR: Failed to allocate abundance_out for species ', nspecies + species_added
            RETURN
          END IF
          s(nspecies + species_added)%abundance_out = 0.0_wp
          
          WRITE(*,'(A,I0,A,A10,A,A10)') 'Added suprathermal species: ', nspecies + species_added, &
            ': ', s(nspecies + species_added)%name, ' based on ', base_name
        END IF
      END IF
    END IF
  END DO
  
  ! Update the species count
  nspecies = nspecies + species_added
  WRITE(*,'(A,I0,A,I0)') 'Added ', species_added, ' suprathermal bulk species, total now: ', nspecies
  
  ! Now print the complete species list after adding suprathermal bulk species
  WRITE(*,'(A,I0)') 'Total species after adding suprathermal bulk species: ', nspecies
  DO i = 1, nspecies
    WRITE(*,'(A,I0,A,A10)') 'Species ', i, ': ', TRIM(s(i)%name)
  END DO
  
  ! Process incoming reaction files that need these species 
  Nlines = 0
  CALL add_reaction_file('radiolysis.dat', 'radiolysis_reactions.out', Nlines, output_unit)
  WRITE(*,'(A,I0,A)') 'Added ', Nlines, ' reactions from radiolysis.dat'
  
  CALL add_reaction_file('class_2_suprathermal.dat', 'class2_reactions.out', Nlines, output_unit)
  WRITE(*,'(A,I0,A)') 'Added ', Nlines, ' reactions from class_2_suprathermal.dat'
  
  CALL add_reaction_file('photo_processes.dat', 'photochemistry_reactions.out', Nlines, output_unit)
  WRITE(*,'(A,I0,A)') 'Added ', Nlines, ' reactions from photo_processes.dat'
  
  ! Handle regular reactions
  Nsup_g = 0
  DO i = 1, nreactions
    IF (r(i)%rtype == 13 .OR. r(i)%rtype == 14) THEN
      ! FIX: Check if r1 and r2 are valid before comparing them
      IF (LEN_TRIM(r(i)%r1) > 0 .AND. LEN_TRIM(r(i)%r2) > 0) THEN
        Nsup_g = Nsup_g + MERGE(1, 2, r(i)%r1 == r(i)%r2)
      ELSE
        Nsup_g = Nsup_g + 2  ! Default to 2 if we can't compare
      END IF
    END IF
  END DO
  
  CALL resize_reaction_array(Nsup_g)
  s_r_counter = 0
  first_suprathermal_react = nreactions + 1
  DO i = first_surfreact, nreactions
    IF (r(i)%rtype == 13 .OR. r(i)%rtype == 14) THEN
      CALL add_suprathermal_reaction(i, s_r_counter)
    END IF
  END DO
  nreactions = nreactions + s_r_counter  ! FIX: Use s_r_counter instead of Nsup_g to be safe
  
  ! FIX: Make sure file_unit is initialized before passing to write_reactions_to_file
  IF (output_unit <= 0) THEN
    CALL find_free_unit(output_unit, MIN_UNIT, MAX_UNIT)
    IF (output_unit == -1) output_unit = MIN_UNIT + 1
  END IF
  CALL write_reactions_to_file(output_unit, 'suprathermal_reactions.out', first_suprathermal_react, nreactions)
  
  ! FIX: Initialize file_unit before passing to add_quenching_reactions
  file_unit = output_unit + 1
  IF (file_unit > MAX_UNIT) file_unit = MIN_UNIT
  CALL find_free_unit(file_unit, MIN_UNIT, MAX_UNIT)
  IF (file_unit == -1) file_unit = MIN_UNIT + 2
  
  ! Now add quenching reactions AFTER we have all species
  CALL add_quenching_reactions(file_unit, output_unit)
  
  ! Final verification
  CALL verify_reactions()
END SUBROUTINE add_suprathermal_bulk_and_reactions

SUBROUTINE write_reactions_to_file(unit_num, filename, start_idx, end_idx)
  INTEGER, INTENT(OUT) :: unit_num
  INTEGER, INTENT(IN) :: start_idx, end_idx
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER :: i, io_stat
  
  ! FIX: Make sure unit_num is valid
  IF (unit_num <= 0 .OR. unit_num >= MAX_UNIT) THEN
    CALL find_free_unit(unit_num, MIN_UNIT, MAX_UNIT)
    IF (unit_num == -1) unit_num = MIN_UNIT + 3  ! Use a default if no units available
  END IF
  
  OPEN(unit_num, FILE=filename, STATUS='REPLACE', IOSTAT=io_stat)
  IF (io_stat == 0) THEN
    DO i = start_idx, end_idx
      ! FIX: Ensure we are within bounds
      IF (i > 0 .AND. i <= SIZE(r)) THEN
        ! FIX: Use the exact same FORMAT statement as used elsewhere in the code (FORMAT 1001)
        WRITE(unit_num, '(1X,I4,1X,I2,1X,2(A10),10X,5(A10),1pE14.5,1X,I1,1X,1pE14.5)', IOSTAT=io_stat) &
          r(i)%idx, r(i)%rtype, &
          r(i)%r1, r(i)%r2, r(i)%p1, &
          r(i)%p2, r(i)%p3, r(i)%p4, &
          r(i)%p5, r(i)%exothermicity, &
          r(i)%exothermicity_known, r(i)%alpha
          
        IF (io_stat /= 0) THEN
          WRITE(*,'(A,I0,A,I0)') 'Warning: Error writing reaction ', i, ' to file: ', io_stat
        END IF
      ELSE
        WRITE(*,'(A,I0,A,I0)') 'Warning: Reaction index ', i, ' out of bounds. Array size: ', SIZE(r)
      END IF
    END DO
    CLOSE(unit_num)
  ELSE
    WRITE(*,'(A,A,A,I0)') 'Error opening file ', TRIM(filename), ' for writing: ', io_stat
  END IF
END SUBROUTINE write_reactions_to_file

END MODULE read_rate06