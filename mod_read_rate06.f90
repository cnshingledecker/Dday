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
  INTEGER :: j, k  ! Added declarations for loop variables
  REAL(KIND=wp) :: a, b, c, apriori_nml
  CHARACTER(LEN=10) :: groundstate, s_name, r1, r2, p1, p2, p3, p4, p5
  CHARACTER(LEN=MAX_FILE_PATH_LEN) :: filename
  TYPE(reaction), DIMENSION(:), ALLOCATABLE :: rtemp
  LOGICAL :: file_exists, is_suprathermal_on
  CHARACTER(LEN=120) :: line
  CHARACTER(LEN=10) :: temp_fields(5)  ! Temporary array for products

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

  ! Allocate arrays and add suprathermal species
  is_suprathermal_on = (suprathermal == 1)
  WRITE(*,'(A,I0)') 'Allocating arrays with nspecies=', nspecies
  CALL allocate_arrays(is_suprathermal_on, max_monolayers)
  IF (is_suprathermal_on) THEN
    WRITE(*,'(A)') 'Adding suprathermal species before main read'
    CALL add_suprathermal_surface_species()
  END IF

  ! Main reading phase
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

  ! Read enthalpies
  CALL read_enthalpias()

  ! Read reactions with dynamic parsing
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
    ! Determine number of fields
    num_fields = 0
    temp_fields = ' '  ! Clear temp array
    READ(line, *, IOSTAT=io_stat) r(ii)%idx, r(ii)%r1, r(ii)%r2, &
      (temp_fields(j), j=1,5), r(ii)%alpha, r(ii)%beta, r(ii)%gamma
    IF (io_stat == 0) THEN
      num_fields = 5  ! 2 products + 3 reals
    ELSE IF (io_stat == -1) THEN
      DO j = 1, 5
        READ(line, *, IOSTAT=io_stat) r(ii)%idx, r(ii)%r1, r(ii)%r2, &
          (temp_fields(k), k=1,j), r(ii)%alpha, r(ii)%beta, r(ii)%gamma
        IF (io_stat == 0) THEN
          num_fields = j  ! Number of products (1 or 2)
          EXIT
        ELSE IF (io_stat /= -1) THEN
          WRITE(*,'(A,I0,A,I0)') 'Failed field count for reaction ', i, ' with IOSTAT=', io_stat
          CALL handle_io_error('counting fields', i)
        END IF
      END DO
    ELSE
      WRITE(*,'(A,I0,A,I0)') 'Failed to parse reaction ', i, ' with IOSTAT=', io_stat
      CALL handle_io_error('reading reaction', i)
    END IF
    ! Assign products based on field count
    IF (num_fields >= 1) THEN
      r(ii)%p1 = temp_fields(1)
      IF (num_fields >= 2) THEN
        r(ii)%p2 = temp_fields(2)
      ELSE
        r(ii)%p2 = ' '
      END IF
      r(ii)%p3 = ' '
      r(ii)%p4 = ' '
      r(ii)%p5 = ' '
    ELSE
      WRITE(*,'(A,I0)') 'Error: No products found for reaction ', i
      CALL handle_io_error('no products', i)
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
  CALL add_bulk_reactions(output_unit)
  WRITE(*,'(A,L1)') 'Added bulk reactions, is_suprathermal_on=', is_suprathermal_on
  IF (is_suprathermal_on) CALL add_suprathermal_bulk_and_reactions(Nsup_g, first_suprathermal_react, output_unit)
  WRITE(*,'(A)') 'Calling write_species_output'
  CALL write_species_output(file_unit)
  WRITE(*,'(A)') 'Finished write_species_output'

  RETURN
1000 FORMAT(1X,I4,1X,2(A10),10X,5(A10),E12.4,1X,F8.2,1X,F8.1)  ! For reference
1001 FORMAT(1X,I4,1X,I2,1X,2(A10),10X,5(A10),1pE14.5,1X,I1,1X,1pE14.5)
END SUBROUTINE read_rate06database

SUBROUTINE assign_species_indices(ii)
  INTEGER, INTENT(IN) :: ii
  INTEGER :: ir1, ir2, ip1, ip2, ip3, ip4, ip5
  ir1 = species_idx(r(ii)%r1)
  WRITE(*,'(A,A10,A,I0)') 'Species ', r(ii)%r1, ' index=', ir1
  ir2 = species_idx(r(ii)%r2)
  WRITE(*,'(A,A10,A,I0)') 'Species ', r(ii)%r2, ' index=', ir2
  ip1 = species_idx(r(ii)%p1)
  WRITE(*,'(A,A10,A,I0)') 'Species ', r(ii)%p1, ' index=', ip1
  ip2 = species_idx(r(ii)%p2)
  WRITE(*,'(A,A10,A,I0)') 'Species ', r(ii)%p2, ' index=', ip2
  ip3 = species_idx(r(ii)%p3)
  WRITE(*,'(A,A10,A,I0)') 'Species ', r(ii)%p3, ' index=', ip3
  ip4 = species_idx(r(ii)%p4)
  WRITE(*,'(A,A10,A,I0)') 'Species ', r(ii)%p4, ' index=', ip4
  ip5 = species_idx(r(ii)%p5)
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

  SUBROUTINE add_suprathermal_bulk_and_reactions(Nsup_g, first_suprathermal_react, output_unit)
    INTEGER, INTENT(OUT) :: Nsup_g, first_suprathermal_react
    INTEGER, INTENT(INOUT) :: output_unit
    INTEGER :: i, ii, s_r_counter, Nlines, io_stat, file_unit
    TYPE(reaction), DIMENSION(:), ALLOCATABLE :: rtemp
    DO i = 1, n_surf_spec
      s(nspecies + i)%name = 'b' // TRIM(s(first_surf_spec + i - 1)%name(2:)) // '*'
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
    Nsup_g = 0
    DO i = 1, nreactions
      IF (r(i)%rtype == 13 .OR. r(i)%rtype == 14) THEN
        Nsup_g = Nsup_g + MERGE(1, 2, r(i)%r1 == r(i)%r2)
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
    nreactions = nreactions + Nsup_g
    CALL write_reactions_to_file(output_unit, 'suprathermal_reactions.out', first_suprathermal_react, nreactions)
    CALL add_reaction_file('radiolysis.dat', 'radiolysis_reactions.out', Nlines, output_unit)
    CALL add_reaction_file('class_2_suprathermal.dat', 'class2_reactions.out', Nlines, output_unit)
    CALL add_quenching_reactions(file_unit, output_unit)
    CALL add_reaction_file('photo_processes.dat', 'photochemistry_reactions.out', Nlines, output_unit)
    CALL verify_reactions()
  END SUBROUTINE add_suprathermal_bulk_and_reactions

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

  SUBROUTINE add_quenching_reactions(file_unit, output_unit)
    INTEGER, INTENT(INOUT) :: file_unit, output_unit
    INTEGER :: ii, jj, io_stat, first_suprathermal_species
    CHARACTER(LEN=10) :: groundstate
    TYPE(reaction), DIMENSION(:), ALLOCATABLE :: rtemp
    INTEGER :: quench_count
    first_suprathermal_species = nspecies - n_surf_spec + 1
    quench_count = nspecies - first_suprathermal_species - 1
    IF (quench_count <= 0) RETURN
    ALLOCATE(rtemp(nreactions + quench_count))
    rtemp(1:SIZE(r)) = r
    CALL find_free_unit(file_unit, MIN_UNIT, MAX_UNIT)
    OPEN(file_unit, FILE='quenching.out', STATUS='REPLACE', IOSTAT=io_stat)
    jj = first_suprathermal_species
    DO ii = nreactions + 1, nreactions + quench_count
      groundstate = s(jj)%name(1:LEN_TRIM(s(jj)%name) - 1)
      rtemp(ii)%r1 = s(jj)%name
      rtemp(ii)%r2 = 'QUENCH'
      rtemp(ii)%p1 = groundstate
      rtemp(ii)%p2 = ' '
      rtemp(ii)%p3 = ' '
      rtemp(ii)%p4 = ' '
      rtemp(ii)%p5 = ' '
      rtemp(ii)%idx = ii
      rtemp(ii)%ir1 = species_idx(rtemp(ii)%r1)
      rtemp(ii)%ir2 = 0
      rtemp(ii)%ip1 = species_idx(groundstate)
      rtemp(ii)%ip2 = 0
      rtemp(ii)%ip3 = 0
      rtemp(ii)%ip4 = 0
      rtemp(ii)%ip5 = 0
      rtemp(ii)%alpha = 1.0_wp
      rtemp(ii)%beta = 1.0_wp
      rtemp(ii)%gamma = 1.0_wp
      rtemp(ii)%rtype = get_rtype(rtemp(ii)%r1, rtemp(ii)%r2)
      rtemp(ii)%exothermicity = 0.0_wp
      rtemp(ii)%exothermicity_known = 0
      WRITE(file_unit, 1000) rtemp(ii)%idx, rtemp(ii)%r1, rtemp(ii)%r2, rtemp(ii)%p1, &
        rtemp(ii)%p2, rtemp(ii)%p3, rtemp(ii)%p4, rtemp(ii)%p5, rtemp(ii)%alpha, &
        rtemp(ii)%beta, rtemp(ii)%gamma
      jj = jj + 1
    END DO
    CLOSE(file_unit)
    CALL resize_and_copy_reactions(rtemp)
    nreactions = nreactions + quench_count
1000 FORMAT(1X,I4,1X,2(A10),10X,5(A10),E12.4,1X,F8.2,1X,F8.1)
  END SUBROUTINE add_quenching_reactions

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
    reac%ir1 = species_idx(reac%r1)
    reac%ir2 = species_idx(reac%r2)
    reac%ip1 = species_idx(reac%p1)
    reac%ip2 = species_idx(reac%p2)
    reac%ip3 = species_idx(reac%p3)
    reac%ip4 = species_idx(reac%p4)
    reac%ip5 = species_idx(reac%p5)
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

  SUBROUTINE write_reactions_to_file(unit_num, filename, start_idx, end_idx)
    INTEGER, INTENT(OUT) :: unit_num
    INTEGER, INTENT(IN) :: start_idx, end_idx
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: i, io_stat
    CALL find_free_unit(unit_num, MIN_UNIT, MAX_UNIT)
    OPEN(unit_num, FILE=filename, STATUS='REPLACE', IOSTAT=io_stat)
    IF (io_stat == 0) THEN
      DO i = start_idx, end_idx
        WRITE(unit_num, 1001) r(i)%idx, r(i)%rtype, r(i)%r1, r(i)%r2, r(i)%p1, &
          r(i)%p2, r(i)%p3, r(i)%p4, r(i)%p5, r(i)%exothermicity, &
          r(i)%exothermicity_known, r(i)%alpha
      END DO
      CLOSE(unit_num)
    END IF
1001 FORMAT(1X,I4,1X,I2,1X,2(A10),10X,5(A10),1pE14.5,1X,I1,1X,1pE14.5)
  END SUBROUTINE write_reactions_to_file

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

END MODULE read_rate06