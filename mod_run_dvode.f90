MODULE run_dvode
USE global_variables
USE global_functions
USE calculate_rates
USE DVODE_F90_M
USE save_results
IMPLICIT NONE
CONTAINS

SUBROUTINE run_dvode_solver
IMPLICIT NONE
INTEGER                             :: i, j, istate, itask,n
REAL(wp), DIMENSION(:), ALLOCATABLE :: y
REAL(wp)                            :: tbegin, tfinish, tt, tcur, tnext
REAL(wp)                            :: tot_surf_ab,tot_bulk_ab
TYPE (VODE_OPTS)                    :: OPTIONS
REAL(wp)                            :: total_abundance, wrt
INTEGER                             :: err_status
CHARACTER(256)                      :: err_iomsg

write(41,*)'time, dtran'
write(42,*)'time, absolute coverage, alpha (fractional coverage), absolute coverage as sum(y(first_surf_spec:nspecies)), bulk abs. cov. by derivative, bult abs. cov. as sum of abundances'

nml = 0
nml_max = 0
s(1:nspecies)%abundance = 0.0d0

DO i =1, init_non_zero
  IF ( MODEL_EXPERIMENT .EQ. 0 ) THEN
    s(species_idx(s_init(i)%name))%abundance = s_init(i)%abundance*gdens/ddens
  ELSE
    IF ( (s(species_idx(s_init(i)%name))%name(1:1) == 'g') .OR. (s(species_idx(s_init(i)%name))%name(1:1) == 'b') ) THEN
!      s(species_idx(s_init(i)%name))%abundance = s_init(i)%abundance*RHO_ICE*(1.0e-8*ICE_THICK) ! equivalent to an ice 1μm x 1μm x ICE_THICK/1e-4 μm
      s(species_idx(s_init(i)%name))%abundance = s_init(i)%abundance*RHO_ICE*(ICE_AREA*ICE_THICK) ! Here we have arbitrarily scaled down the number of ice molecules to a 1e-14 cm2 * IceTHick volume 
!      s(species_idx(s_init(i)%name))%abundance = s_init(i)%abundance*RHO_ICE*1.0e-12 ! 1 um^3 in cm^3
    ElSE
      s(species_idx(s_init(i)%name))%abundance = s_init(i)%abundance*gdens/ddens
    ENDIF
  ENDIF

  IF (s(species_idx(s_init(i)%name))%name == 'G0'.OR.s(species_idx(s_init(i)%name))%name == 'G-'.OR.s(species_idx(s_init(i)%name))%name == 'G+') s(species_idx(s_init(i)%name))%abundance = s(species_idx(s_init(i)%name))%abundance*ddens/gdens
!  print*,s(species_idx(s_init(i)%name))%name,s(species_idx(s_init(i)%name))%abundance
  IF (s(species_idx(s_init(i)%name))%name(1:1) == 'g') nml = nml + s(species_idx(s_init(i)%name))%abundance
  IF (s(species_idx(s_init(i)%name))%name(1:1) == 'b') nml = nml + s(species_idx(s_init(i)%name))%abundance
  ! Question: why not also add 'b' species? - CNS
  IF (s(species_idx(s_init(i)%name))%name(1:4) == 'bH2O') initial_water = s(species_idx(s_init(i)%name))%abundance
  IF (s(species_idx(s_init(i)%name))%name(1:3) == 'bO2') initial_oxygen = s(species_idx(s_init(i)%name))%abundance
ENDDO

nml = aint(nml/nsites)
nml_old = -1

itask = 1
istate = 1
tbegin = tstart*year
tfinish   = tend*year
tt = 0
tcur = tbegin

! y contains the abundances of all species, with 1 entry per species
ALLOCATE(y(nspecies))
y(:) = 0.0d0
y(1:nspecies) = s(1:nspecies)%abundance
PRINT *, "Sum(Y) =",sum(y(1:nspecies))
PRINT *, "Initial number of monolayers = ",nml

total_atoms = 0
DO i = first_surf_spec,nspecies
  IF ( TRIM(s(i)%name) .NE. "e-" .AND. TRIM(s(i)%name) .NE. "ge-" .AND. TRIM(s(i)%name).NE."be-" ) THEN
    WRITE (*,'(A,A,ES10.4,A,ES10.4)') s(i)%name, "has an s-array abundance of ", s(i)%abundance, " and a y-array abundance of ",y(i)
    total_atoms = total_atoms + s(i)%abundance*s(i)%natoms
  ENDIF
ENDDO

WRITE (*,'(A,ES10.4,A)') "There are a total of ", total_atoms, " atoms at the beginning"


OPTIONS = SET_OPTS(MXSTEP=500000, RELERR=RTOL, ABSERR=ATOL, METHOD_FLAG=22)

CALL open_analytics_files

DO i = 1, timesteps

  tnext = tcur*10.d0**(log10(tend/tstart)/timesteps)

!  DO n = 1,nspecies
!    IF ( ISNAN(y(n)) ) THEN
!      PRINT *, s(n)%name, " is NaN"
!      y(n) = 0.0d0
!      s(n)%abundance = 0.0d0
!    ENDIF
!  ENDDO

  ! Do chemistry this timestep
  CALL DVODE_F90(re,nspecies,Y,tt,tnext,ITASK,ISTATE,OPTIONS)



  ! Save results
  CALL save_analytics(nspecies, Y, tt, i)

  ! Save time
  timesteps_out(i) = tnext

  ! Now update the species abundances in the struct
  DO j = 1, nspecies
!    PRINT *, "'",s(j)%name,"' has length=",LEN_TRIM(s(j)%name)
!    IF ( s(j)%name(1:1) .EQ. '' ) THEN
!      PRINT *, "This name is empty!!"
!      CALL EXIT()
!    ENDIF
    s(j)%abundance_out(i) = y(j)
!    PRINT *, "****************************************"
  ENDDO

  IF ( MODEL_EXPERIMENT .EQ. 1 ) THEN
    IF ( tcur .GT. TEND*3.155e7) GOTO 93
  ENDIF
  tcur = tnext

!  CALL save_results_bulk

ENDDO

93 CONTINUE

DO j = 1, nspecies
  IF ( ( ( s(j)%name(1:1) .EQ. 'b' ) .OR. ( s(j)%name(1:1) .EQ. 'b' ) ) .AND. &
    ( s(j)%frac_abundance .GT. 1.0e-30 ) ) THEN
    PRINT *, s(j)%name, s(j)%frac_abundance
  ENDIF
ENDDO

END SUBROUTINE run_dvode_solver

SUBROUTINE run_dvode_solver_sparse
IMPLICIT NONE

INTEGER                             :: i, j, istate, itask, nia, nja, col, row
REAL(wp), DIMENSION(:), ALLOCATABLE :: y
INTEGER, DIMENSION(:), ALLOCATABLE  :: ia, ja
REAL(wp)                            :: tbegin, tfinish, tt, tcur, tnext
TYPE (VODE_OPTS)                    :: OPTIONS

write(41,*)'time, dtran'

DO i =1, init_non_zero
  s(species_idx(s_init(i)%name))%abundance = s_init(i)%abundance*gdens/ddens
  IF (s(species_idx(s_init(i)%name))%name == 'G0'.OR.s(species_idx(s_init(i)%name))%name == 'G-'.OR.s(species_idx(s_init(i)%name))%name == 'G+') s(species_idx(s_init(i)%name))%abundance = s(species_idx(s_init(i)%name))%abundance*ddens/gdens
  print*,s(species_idx(s_init(i)%name))%name,s(species_idx(s_init(i)%name))%abundance
  IF (s(species_idx(s_init(i)%name))%name(1:1) == 'g') nml = nml + s(species_idx(s_init(i)%name))%abundance
ENDDO

nml = aint(nml/nsites)
nml_old = -1

itask = 1
istate = 1
tbegin = tstart*year
tfinish   = tend*year
tt = 0.d0
tcur = tbegin

ALLOCATE(y(nspecies))
y(1:nspecies) = s(1:nspecies)%abundance

OPTIONS = SET_OPTS(MXSTEP=500000, RELERR=RTOL, ABSERR=ATOL, METHOD_FLAG=126, MA28_RPS=.TRUE.)

CALL open_analytics_files

DO i = 1, timesteps

  tnext = tcur*10.d0**(log10(tend/tstart)/timesteps)

  CALL save_analytics(nspecies, Y, tt, i)

  timesteps_out(i) = tnext
  DO j = 1, nspecies
    s(j)%abundance_out(i) = y(j)
  ENDDO
  tcur = tnext
!    if (tcur>3.155d7*1.0d6) CALL save_results_semenov
!    CALL save_results_bulk

ENDDO

CALL close_analytics_files

END SUBROUTINE run_dvode_solver_sparse

SUBROUTINE run_dvode_solver_sparse_mre
IMPLICIT NONE

INTEGER                             :: i, j, istate, itask, nia, nja, col, row
REAL(wp), DIMENSION(:), ALLOCATABLE :: y
INTEGER, DIMENSION(:), ALLOCATABLE  :: ia, ja
REAL(wp)                            :: tbegin, tfinish, tt, tcur, tnext
TYPE (VODE_OPTS)                    :: OPTIONS

DO i =1, init_non_zero
  s(species_idx(s_init(i)%name))%abundance = s_init(i)%abundance*gdens/ddens
  IF (s(species_idx(s_init(i)%name))%name == 'G0'.OR.s(species_idx(s_init(i)%name))%name == 'G-'.OR.s(species_idx(s_init(i)%name))%name == 'G+') s(species_idx(s_init(i)%name))%abundance = s(species_idx(s_init(i)%name))%abundance*ddens/gdens
ENDDO

itask = 1
istate = 1
tbegin = tstart*year
tfinish   = tend*year
tt = 0
tcur = tbegin

ALLOCATE(y(nspecies))
y(1:nspecies) = s(1:nspecies)%abundance

OPTIONS = SET_OPTS(MXSTEP=500000, RELERR=RTOL, ABSERR=ATOL, METHOD_FLAG=22, MA28_RPS=.TRUE.)

CALL open_analytics_files

DO i = 1, timesteps

  tnext = tcur*10.d0**(log10(tend/tstart)/timesteps)

!  CALL DVODE_F90(mre,nspecies,Y,tt,tnext,ITASK,ISTATE,OPTIONS)
  CALL save_analytics(nspecies, Y, tt, i)

  timesteps_out(i) = tnext
  DO j = 1, nspecies
    s(j)%abundance_out(i) = y(j)
  ENDDO
  tcur = tnext

ENDDO

CALL close_analytics_files

END SUBROUTINE run_dvode_solver_sparse_mre

SUBROUTINE re(neq, t, y, ydot)
INTEGER, INTENT (IN)     :: neq
INTEGER                  :: i, j, k, nml_int,bulk_idx
REAL(wp)                 :: rr, transit, cov, diff_s2m, diff_m2s, diff_m2s_tot, tot_diff, dtran_fd, dtran_sr, diff_s2m_tot, dtran_tmp, d3, d4, dtran_chb
REAL(wp)                 :: dtran_freeze,dtran_desorb,dtran_chemdes
REAL(wp)                 :: cd_p1, cd_p2, cd_p3, cd_p4, cd_p5, effmass, non_h2o_surf_frac
REAL(wp)                 :: total_abundance, wrt
REAL(wp), INTENT (IN) :: t
REAL(wp), DIMENSION(neq), INTENT(IN) :: y
REAL(wp), DIMENSION(neq), INTENT(OUT) :: ydot
character*10 :: sss
logical :: file_exists
REAL(wp) :: f1,f2,tot_surf_ab,tot_bulk_ab,total_s,temp_element
REAL(wp) :: temp_atoms,temp_atoms_y, nml_s_array

IF (delta_rho==1 .OR. delta_t==1) CALL calc_rates(t)




1001 FORMAT(1X,I4,1X,I2,1X,2(A10),10X,5(A10),1pE14.5,1X,1pE14.5,1x,1pE14.5)
1002 FORMAT(1X,I4,1X,I2,1X,2(A10),10X,5(A10),1pE14.5,1X,1pE14.5,1x,1pE14.5, 1x, 1pE14.5)
1003 FORMAT(1X,I4,1X,I2,1X,2(A10),10X,5(A10),7(1pE14.5,1X))

ydot(:) = 0.0d0

! This gives the number of monolayers
cov  = sum(y(first_surf_spec:first_bulk_spec-1))/nsites !(nsites*n_s_ml)

non_h2o_surf_frac = 0.0d0
IF ( MODEL_EXPERIMENT .EQ. 0 ) THEN
  non_h2o_surf_frac = 1.0d0 - y(species_idx('gH2O      '))/(sum(y(first_surf_spec:first_bulk_spec-1)) + 1.0d-99)
ENDIF

! smol is the sum of the abundance of surface species
smol = sum(y(first_surf_spec:first_bulk_spec-1))


! bmol is the sum of the abundance of bulk species
bmol = sum(y(first_bulk_spec:nspecies))

dtran_fd = 0.0d0
dtran_sr = 0.0d0
dtran_chb = 0.0d0

dtran_freeze = 0.0d0
dtran_desorb = 0.0d0
dtran_chemdes = 0.0d0

! Get the total abundance of bulk species
tot_bulk_ab = 0.0
DO j = 1, nspecies
  IF  ( s(j)%name(1:1) .EQ. 'b' ) THEN
    tot_bulk_ab = tot_bulk_ab + y(j)
  ENDIF
ENDDO


! Get the total abundance of surface species
tot_surf_ab = 0.0
DO j = 1, nspecies
  IF  ( s(j)%name(1:1) .EQ. 'g' ) THEN
    tot_surf_ab = tot_surf_ab + y(j)
  ENDIF
ENDDO

! MODIFICATION: SET smol = tot_surf_ab and bmol = tot_bulk_ab
smol = tot_surf_ab
bmol = tot_bulk_ab


DO j = 1, nreactions
  rr = 0.0d0
  IF (r(j)%ir2==0) THEN
    rr = r(j)%rate*y(r(j)%ir1)
  ELSE IF ( (r(j)%rtype .EQ. 17) .OR. (r(j)%rtype .EQ. 18) .OR. (r(j)%rtype .EQ. 19) .OR. (r(j)%rtype .EQ. 20)) THEN
    rr = r(j)%rate*y(r(j)%ir1)
  ELSE
    rr = r(j)%rate*y(r(j)%ir1)*y(r(j)%ir2)
  ENDIF

    !Bulk reactions:
    IF ( r(j)%rtype==14 .OR. r(j)%rtype==16 ) THEN
        IF (bmol<nsites) THEN
!          IF ( ANY(radical_names .EQ. r(j)%r1 ) .OR. ANY(radical_names .EQ. r(j)%r2 ) ) THEN
!            CONTINUE
!          ELSE
            rr = rr/nsites
!          ENDIF
        ELSE
!          IF ( ANY(radical_names .EQ. r(j)%r1 ) .OR. ANY(radical_names .EQ. r(j)%r2 ) ) THEN
!            CONTINUE
!          ELSE
            rr = rr/bmol
!          ENDIF
        ENDIF
    ENDIF

!    IF ( FAST_BULK .EQ. 1 ) THEN
!      IF (  (r(j)%rtype .EQ. 16) .OR. &
!        ( ANY(radical_names .EQ. r(j)%r1) .OR. ANY(radical_names .EQ. r(j)%r2) ) ) THEN
!        f1 = y(r(j)%ir1)/tot_bulk_ab
!        f2 = y(r(j)%ir2)/tot_bulk_ab
!        rr = rr*(f1*f2)
       !          IF ( ( f1 .GT. 0.0 ) .AND. ( f2 .GT. 0.0 ) ) THEN
        !            PRINT *, r(j)%r1," with fraction of ",f1
        !            PRINT *, r(j)%r2," with fraction of ",f2
        !            PRINT *, "Combined probability = ",f1*f2
        !            PRINT *, "Rate coefficient     = ",r(j)%rate*f1*f2
        !            PRINT *, "Rate                 = ",rr
        !            PRINT *, "********************************"
        !          ENDIF
!      ENDIF
!    ENDIF


!  IF (((s(r(j)%ir1)%name(1:4) .EQ. 'bSO ') .OR.  (s(r(j)%ir2)%name(1:4) .EQ. 'bSO ')) &
!    .AND. (rr .GT. 0.0e0)) THEN
!  IF (r(j)%rtype==17 .AND. ((s(r(j)%ir1)%name(1:4) .EQ. 'bSO '))) THEN
!  IF ( r(j)%rtype==18 .AND. (((s(r(j)%ir1)%name(1:5) .EQ. 'bSO3*')))) THEN
!  IF ((r(j)%rtype==16 .AND. (((s(r(j)%ir1)%name(1:4) .EQ. 'bOH*')) .AND. &
!                            ((s(r(j)%ir2)%name(1:5) .EQ. 'bH2O ')))) .OR. &
!     (r(j)%rtype==16 .AND. (((s(r(j)%ir1)%name(1:4) .EQ. 'bOH '))   .AND. &
!                            ((s(r(j)%ir2)%name(1:5) .EQ. 'bH2O*'))))) THEN
!  IF ((s(r(j)%ir1)%name(1:3) .EQ. 'bO ') .AND.  (s(r(j)%ir2)%name(1:3) .EQ. 'bO ')) THEN
!  IF ( r(j)%rtype .EQ. 18 ) THEN
!    PRINT '(A10,A,A10,A,A10,A,A10,A,A10)', r(j)%r1," + ",r(j)%r2," -> ",r(j)%p1," + ",r(j)%p2," + ",r(j)%p3
!    PRINT *, "idx = ",r(j)%idx
!    PRINT *, "rtype =",r(j)%rtype
!    PRINT '(A,ES10.4)', "Reaction rate = ",rr
!    PRINT '(A,A,A,ES10.4)', "Abundance of ",r(j)%r1," =",y(r(j)%ir1)
!    PRINT '(A,A,A,ES10.4)', "Abundance of ",r(j)%r2," =",y(r(j)%ir2)
!    PRINT *, "exothermicity_known = ",r(j)%exothermicity_known
!    PRINT *, "The reactants/products of this reaction are:"
!    PRINT *, "r1 = ",r(j)%r1,"and ir1 = ",r(j)%ir1
!    PRINT *, "r2 = ",r(j)%r2,"and ir2 = ",r(j)%ir2
!    PRINT *, "p1 = ",r(j)%p1,"and ip1 = ",r(j)%ip1
!    PRINT *, "p2 = ",r(j)%p2,"and ip2 = ",r(j)%ip2
!    PRINT *, "p3 = ",r(j)%p3,"and ip3 = ",r(j)%ip3
!    PRINT *, "p4 = ",r(j)%p4,"and ip4 = ",r(j)%ip4
!    PRINT *, "p5 = ",r(j)%p5,"and ip5 = ",r(j)%ip5
!    PRINT '(A,ES10.4,A,ES10.4,A,ES10.4,A,ES10.4)', "Rate coefficient= ",r(j)%rate," : alpha = ",r(j)%alpha," beta =" ,r(j)%beta," gamma = ",r(j)%gamma
!    PRINT *, "Exothermicity = ", r(j)%exothermicity
!    PRINT *, "***********************************"
    IF ( ISNAN(rr) ) THEN
      PRINT *, "RR = NaN"
      CALL EXIT()
    ENDIF
!
!  ENDIF


!

!  IF (((s(r(j)%ir1)%name(1:4) .EQ. 'bSO ') .OR.  (s(r(j)%ir2)%name(1:4) .EQ. 'bSO ')) &
!    .AND. (rr .GT. 0.0e0)) THEN
!  IF (r(j)%rtype==17 .AND. ((s(r(j)%ir1)%name(1:4) .EQ. 'bSO '))) THEN
!  IF ( r(j)%rtype==18 .AND. (((s(r(j)%ir1)%name(1:5) .EQ. 'bSO3*')))) THEN
!  IF (r(j)%rtype==18 .AND. (((s(r(j)%ir1)%name(1:4) .EQ. 'bSO2*')) .OR. &
!                            ((s(r(j)%ir1)%name(1:4) .EQ. 'bSO3*')) .OR. &
!                            ((s(r(j)%ir1)%name(1:4) .EQ. 'bSO3*')))) THEN
!    PRINT *, r(j)%r1," + ",r(j)%r2," -> ",r(j)%p1," + ",r(j)%p2," + ",r(j)%p3
!    PRINT *, "Reaction rate = ",rr
!    PRINT *, "Rate coefficient= ",r(j)%rate
!    PRINT *, "Abundance of ",r(j)%r1," =",y(r(j)%ir1)/1.0e20
!    PRINT *, "Abundance of ",r(j)%r2," =",y(r(j)%ir2)
!    PRINT *, "***********************************"
!  ENDIF

  ydot(r(j)%ir1) = ydot(r(j)%ir1) - rr
  IF (r(j)%ir2/=0) ydot(r(j)%ir2) = ydot(r(j)%ir2) - rr
  ydot(r(j)%ip1) = ydot(r(j)%ip1) + rr
  IF (r(j)%ip2/=0) ydot(r(j)%ip2) = ydot(r(j)%ip2) + rr
  IF (r(j)%ip3/=0) ydot(r(j)%ip3) = ydot(r(j)%ip3) + rr
  IF (r(j)%ip4/=0) ydot(r(j)%ip4) = ydot(r(j)%ip4) + rr
  IF (r(j)%ip5/=0) ydot(r(j)%ip5) = ydot(r(j)%ip5) + rr


  !REACTIVE DESORPTION ACCORDING TO MINISSALE&DULIEU:
  IF (des_reactive_type == 2) THEN
      IF (((r(j)%rtype == 13) .OR. (r(j)%rtype == 15)) .AND. r(j)%p1(1:1) == 'g') THEN
          cd_p1 = des_reactive
          cd_p2 = des_reactive
          cd_p3 = des_reactive
          cd_p4 = des_reactive
          cd_p5 = des_reactive

          IF (r(j)%exothermicity_known == 1) THEN

              effmass = ((s(r(j)%ip1)%weight-effsurfmass)/(s(r(j)%ip1)%weight+effsurfmass))**2.0d0
              cd_p1 = dexp( -s(r(j)%ip1)%natoms*3.0d0 * s(r(j)%ip1)%edes / r(j)%exothermicity / effmass ) * non_h2o_surf_frac
              IF (cd_p1>1.0d0) cd_p1 = 0.0d0
              ydot(r(j)%ip1) = ydot(r(j)%ip1) - cd_p1*rr
              ydot(s(r(j)%ip1)%gas_idx) = ydot(s(r(j)%ip1)%gas_idx) + cd_p1*rr
              rd_v2_terms(1,j) = cd_p1*rr
              rd_v2_terms(1,j+1) = cd_p1*rr

              IF (r(j)%ip2/=0) THEN
                  effmass = ((s(r(j)%ip2)%weight-effsurfmass)/(s(r(j)%ip2)%weight+effsurfmass))**2.0d0
                  cd_p2 = dexp( -s(r(j)%ip2)%natoms*3.0d0 * s(r(j)%ip2)%edes / r(j)%exothermicity / effmass ) * non_h2o_surf_frac
                  IF (cd_p2>1.0d0) cd_p2 = 0.0d0
                  ydot(r(j)%ip2) = ydot(r(j)%ip2) - cd_p2*rr
                  ydot(s(r(j)%ip2)%gas_idx) = ydot(s(r(j)%ip2)%gas_idx) + cd_p2*rr
                  rd_v2_terms(2,j) = cd_p2*rr
                  rd_v2_terms(2,j+1) = cd_p2*rr
              ENDIF

              IF (r(j)%ip3/=0) THEN
                  effmass = ((s(r(j)%ip3)%weight-effsurfmass)/(s(r(j)%ip3)%weight+effsurfmass))**2.0d0
                  cd_p3 = dexp( -s(r(j)%ip3)%natoms*3.0d0 * s(r(j)%ip3)%edes / r(j)%exothermicity / effmass ) * non_h2o_surf_frac
                  IF (cd_p3>1.0d0) cd_p3 = 0.0d0
                  ydot(r(j)%ip3) = ydot(r(j)%ip3) - cd_p3*rr
                  ydot(s(r(j)%ip3)%gas_idx) = ydot(s(r(j)%ip3)%gas_idx) + cd_p3*rr
                  rd_v2_terms(3,j) = cd_p3*rr
                  rd_v2_terms(3,j+1) = cd_p3*rr
              ENDIF

              IF (r(j)%ip4/=0) THEN
                  effmass = ((s(r(j)%ip4)%weight-effsurfmass)/(s(r(j)%ip4)%weight+effsurfmass))**2.0d0
                  cd_p4 = dexp( -s(r(j)%ip4)%natoms*3.0d0 * s(r(j)%ip4)%edes / r(j)%exothermicity / effmass ) * non_h2o_surf_frac
                  IF (cd_p4>1.0d0) cd_p4 = 0.0d0
                  ydot(r(j)%ip4) = ydot(r(j)%ip4) - cd_p4*rr
                  ydot(s(r(j)%ip4)%gas_idx) = ydot(s(r(j)%ip4)%gas_idx) + cd_p4*rr
                  rd_v2_terms(4,j) = cd_p4*rr
                  rd_v2_terms(4,j+1) = cd_p4*rr
              ENDIF

              IF (r(j)%ip5/=0) THEN
                  effmass = ((s(r(j)%ip5)%weight-effsurfmass)/(s(r(j)%ip5)%weight+effsurfmass))**2.0d0
                  cd_p5 = dexp( -s(r(j)%ip5)%natoms*3.0d0 * s(r(j)%ip5)%edes / r(j)%exothermicity / effmass ) * non_h2o_surf_frac
                  IF (cd_p5>1.0d0) cd_p5 = 0.0d0
                  ydot(r(j)%ip5) = ydot(r(j)%ip5) - cd_p5*rr
                  ydot(s(r(j)%ip5)%gas_idx) = ydot(s(r(j)%ip5)%gas_idx) + cd_p5*rr
                  rd_v2_terms(5,j) = cd_p5*rr
                  rd_v2_terms(5,j+1) = cd_p5*rr
              ENDIF

            ELSE
              ydot(r(j)%ip1) = ydot(r(j)%ip1) - cd_p1*rr
              ydot(s(r(j)%ip1)%gas_idx) = ydot(s(r(j)%ip1)%gas_idx) + cd_p1*rr
                rd_v2_terms(1,j) = cd_p1*rr
                rd_v2_terms(1,j+1) = cd_p1*rr

              IF (r(j)%ip2/=0) ydot(r(j)%ip2) = ydot(r(j)%ip2) - cd_p2*rr
              IF (r(j)%ip2/=0) ydot(s(r(j)%ip2)%gas_idx) = ydot(s(r(j)%ip2)%gas_idx) + cd_p2*rr
                rd_v2_terms(2,j) = cd_p2*rr
                rd_v2_terms(2,j+1) = cd_p2*rr

              IF (r(j)%ip3/=0) ydot(r(j)%ip3) = ydot(r(j)%ip3) - cd_p3*rr
              IF (r(j)%ip3/=0) ydot(s(r(j)%ip3)%gas_idx) = ydot(s(r(j)%ip3)%gas_idx) + cd_p3*rr
                rd_v2_terms(3,j) = cd_p3*rr
                rd_v2_terms(3,j+1) = cd_p3*rr

              IF (r(j)%ip4/=0) ydot(r(j)%ip4) = ydot(r(j)%ip4) - cd_p4*rr
              IF (r(j)%ip4/=0) ydot(s(r(j)%ip4)%gas_idx) = ydot(s(r(j)%ip4)%gas_idx) + cd_p4*rr
                rd_v2_terms(4,j) = cd_p4*rr
                rd_v2_terms(4,j+1) = cd_p4*rr

              IF (r(j)%ip5/=0) ydot(r(j)%ip5) = ydot(r(j)%ip5) - cd_p5*rr
              IF (r(j)%ip5/=0) ydot(s(r(j)%ip5)%gas_idx) = ydot(s(r(j)%ip5)%gas_idx) + cd_p5*rr
                rd_v2_terms(5,j) = cd_p5*rr
                rd_v2_terms(5,j+1) = cd_p5*rr
          ENDIF
          if (rde==0) write(myunit,'(7a10, 6x)')r(j)%r1, r(j)%r2, r(j)%p1, r(j)%p2, r(j)%p3, r(j)%p4, r(j)%p5
          if (rde==0) write(myunit,'(5(1pe16.8))') cd_p1, cd_p2, cd_p3, cd_p4, cd_p5
      ENDIF

      if (rde==0 .and. j==nreactions) then
          rde=1
          close(myunit)
      endif

  ENDIF

  !DIAGNOSTIC TERMS!

    ! Collect freeze rate
    IF (r(j)%rtype==11) THEN
      dtran_fd = dtran_fd + rr
      dtran_freeze = dtran_freeze + rr
    ENDIF

    ! Collect desorb rate
    IF (r(j)%rtype==12) THEN
      dtran_fd = dtran_fd - rr
      dtran_desorb = dtran_desorb + rr
    ENDIF

    ! Collect chemi-desorption rate
    IF ((r(j)%rtype==13 .or. r(j)%rtype==15).AND.(r(j)%p1(1:1).NE.'g')) THEN
        dtran_chemdes = dtran_chemdes + rr
    ENDIF

    if ( &
        r(j)%rtype==13 .or. &
        r(j)%rtype==15 .OR. &
        (r(j)%rtype==2 .and. r(j)%r1(1:1)=='g') .or. &
        (r(j)%rtype==3 .and.  r(j)%r1(1:1)=='g') .OR. &
        (r(j)%rtype==17 .and.  r(j)%r1(1:1)=='g') .OR. &
        (r(j)%rtype==18 .and.  r(j)%r1(1:1)=='g') .OR. &
        (r(j)%rtype==21 .and.  r(j)%r1(1:1)=='g') &
        ) then
        if (r(j)%r1(1:1)=='g') dtran_sr = dtran_sr - rr
        if (r(j)%r2(1:1)=='g') dtran_sr = dtran_sr - rr
        if (r(j)%p1(1:1)=='g') dtran_sr = dtran_sr + rr
        if (r(j)%p2(1:1)=='g') dtran_sr = dtran_sr + rr
        if (r(j)%p3(1:1)=='g') dtran_sr = dtran_sr + rr
        if (r(j)%p4(1:1)=='g') dtran_sr = dtran_sr + rr
        if (r(j)%p5(1:1)=='g') dtran_sr = dtran_sr + rr
    endif

  !END OF DIAGNOSTIC TERMS!

    if ((&
         r(j)%r1(1:1)=='g' .or. &
         r(j)%r2(1:1)=='g' .or. &
         r(j)%p1(1:1)=='g' .or. &
         r(j)%p2(1:1)=='g' .or. &
         r(j)%p3(1:1)=='g' .or. &
         r(j)%p4(1:1)=='g' .or. &
         r(j)%p5(1:1)=='g') .and. &
         (&
         r(j)%rtype/=2 .and. &
         r(j)%rtype/=3 .and. &
         r(j)%rtype/=11 .and. &
         r(j)%rtype/=12 .and. &
         r(j)%rtype/=13 .AND. &
         r(j)%rtype/=15 .AND. &
         r(j)%rtype/=16 .AND. &
         r(j)%rtype/=17 .AND. &
         r(j)%rtype/=18 .AND. &
         r(j)%rtype/=19 .AND. &
         r(j)%rtype/=21 .AND. &
         r(j)%rtype/=20       &
         ) ) pause

ENDDO




!Calculating dtran - total transition rate between bulk and surface due to chemical processes
dtran = 0.0d0 !dtran === (dn_s/dt)0

! This should be over only the surface species? - CNS
!DO j = first_surf_spec, first_surf_spec + n_surf_spec - 1
IF ( s(first_bulk_spec)%name(1:1) .NE. 'b' ) THEN
  PRINT *, s(first_bulk_spec)%name, " is not a bulk species when it should be in mod_run_dvode"
  CALL EXIT()
ENDIF

DO j = first_surf_spec, first_bulk_spec - 1
  ! Loop over all surface species
  IF ( s(j)%name(1:1) .NE. "g" ) THEN
    PRINT *, "Looping over what is not a surface species at line 593"
    CALL EXIT()
  ENDIF
!    PRINT *, s(j)%name," at j=",j
    dtran = dtran + ydot(j)
ENDDO

! This should be from the first bulk species to nspecies? - CNS
!do j = first_surf_spec+n_surf_spec, first_surf_spec+n_surf_spec+n_surf_spec-1
do j = first_bulk_spec, nspecies 
  ! Loop over all bulk species
  IF ( s(j)%name(1:1) .NE. "b" ) THEN
    PRINT *, "Looping over what is not a bulk species at line 605"
    CALL EXIT()
  ENDIF
  dtran_chb = dtran_chb + ydot(j)
enddo

!Building transition terms for species on surface and in bulk
! CNS: Updated loop ranges, which should only be over surface species
transit = 0.0d0
dtran_tmp = 0.0d0
!DO j = first_surf_spec, first_surf_spec + n_surf_spec - 1
DO j = first_surf_spec, first_bulk_spec - 1
  bulk_idx = get_bulk_idx(s(j)%name)
  ! Loop over all surface species
  IF ( s(j)%name(1:1) .NE. "g" ) THEN
    PRINT *, "Looping over what is not a surface species at line 610"
    CALL EXIT()
  ENDIF
  IF (dtran>=0.0d0) THEN
    ! If the change in the number of surface species increases, then it may be
    !necessary to move some to bulk
    IF (smol>0.0d0) transit = dtran*y(j)/(nsites*n_s_ml)
    if (tdust>150.0d0 .AND. cov<5.0d-3) transit = 0.0d0
    ydot(j) = ydot(j) - transit
    ! The old assumption here was that the species at j is the surface
    ! counterpart of the species at j+n_surf_spec (bulk?)
    ! C. N. Shingledecker: Now we assign the bulk_idx using the function in mod_global_functions.f90
    ydot(bulk_idx) = ydot(bulk_idx) + transit
    IF ( s(j)%name(2:LEN_TRIM(s(j)%name)).NE.s(bulk_idx)%name(2:LEN_TRIM(s(j)%name)) ) THEN
      PRINT *, s(j)%name(2:LEN_TRIM(s(j)%name)), " not equal to ",s(bulk_idx)%name(2:LEN_TRIM(s(j)%name))
      CALL EXIT()
    ENDIF
    dtran_tmp = dtran_tmp + transit
  ELSE
    ! If the change in the number of surface species is negative, then it may be
    ! necessary to move some species from the bulk to the surface
    IF (bmol>0.0d0) transit = dtran*dmin1(bmol,smol)/smol*y(j+n_surf_spec)/bmol !y(nspecies+3)
    if (tdust>150.0d0 .AND. cov<5.0d-3) transit = 0.0d0
    ydot(j) = ydot(j) - transit
    ydot(bulk_idx) = ydot(bulk_idx) + transit
    dtran_tmp = dtran_tmp + transit
  ENDIF
ENDDO


! d3 is supposedly the total of the rate of change in bulk species, after
! accounting for mantle growth
! CNS: Updated range of loop, which should only be over bulk species
!d3 = sum(ydot(first_surf_spec+n_surf_spec:first_surf_spec+n_surf_spec+n_surf_spec-1))
d3 = sum(ydot(first_bulk_spec:nspecies))

!Building transition terms due to bulk diffusion (swapping):
diff_m2s_tot = 0.0d0
diff_s2m_tot = 0.0d0

! Loop over all (and only) bulk species
! CNS: Updated range of loop, which should only be over bulk species
!DO i = first_surf_spec + n_surf_spec, first_surf_spec + n_surf_spec + n_surf_spec - 1 !Loop through bulk species
DO i = first_bulk_spec,nspecies !Loop through bulk species
    IF ( s(i)%name(1:1) .NE. "b" ) THEN
      PRINT *, "Looping over what is not a bulk species at line 651"
      CALL EXIT()
    ENDIF
    diff_m2s = 0.0d0
    IF (bmol>0.0d0) diff_m2s = dmin1(1.0d0, smol/bmol)*y(i)*s(i)%bdiffrate * bulk_diff_slowdown
    ydot(i-n_surf_spec) = ydot(i-n_surf_spec) + diff_m2s
    ydot(i) = ydot(i) - diff_m2s
    diff_m2s_tot = diff_m2s_tot + diff_m2s
ENDDO

! Loop over all (and only) bulk species, again
! CNS: Updated range of loop, which should only be over bulk species
!temp_atoms = 0
!temp_atoms_y = 0
!PRINT *, "*******************************************"
!PRINT *, " At line 650... should only be bulk species"
!PRINT *, "*******************************************"
DO i = first_bulk_spec, nspecies !Second loop through bulk species
    IF ( s(i)%name(1:1) .NE. "b" ) THEN
      PRINT *, "Looping over what is not a bulk species at line 670"
      CALL EXIT()
    ENDIF

    ! Sanity check on abundances of bulk species
!    WRITE (*,'(A,I3,A,I3,A,A,A,ES9.2,A,ES9.2,A,F5.1,A)') "At ", i, "of ", nspecies, ": ",s(i)%name, "has an s-array abundance of ", s(i)%abundance, " and a y-array abundance of ",y(i), " and ", ((s(i)%abundance*s(i)%natoms)/total_atoms)*100.0, "% of the initial oxygen"

!    IF ( TRIM(s(i)%name) .NE. "e-" .AND. TRIM(s(i)%name) .NE. "ge-" .AND. TRIM(s(i)%name).NE."be-" ) THEN
!      temp_atoms = temp_atoms + s(i)%abundance*s(i)%natoms
!      temp_atoms_y = temp_atoms_y + y(i)*s(i)%natoms
!    ENDIF

    diff_s2m = 0.0d0
    IF (smol>0.0d0) diff_s2m = y(i-n_surf_spec)/smol*diff_m2s_tot*cov
    ydot(i-n_surf_spec) = ydot(i-n_surf_spec) - diff_s2m
    ydot(i) = ydot(i) + diff_s2m
    diff_s2m_tot = diff_s2m_tot + diff_s2m
ENDDO

!WRITE( *,'(A,ES9.2,A,ES9.2)') "Sum of grain s-array = ",SUM(s(first_surf_spec:nspecies)%abundance), " and sum of grain y-array = ",SUM(y(first_surf_spec:nspecies))
!WRITE( *,'(A,F5.1,A)') "There are ", (temp_atoms/total_atoms)*100.0, "% of the initial atoms left in the s-array"
!WRITE( *,'(A,F5.1,A)') "There are ", (temp_atoms_y/total_atoms)*100.0, "% of the initial atoms left in the y-array"
!PRINT *, "***************************************"



! Changed the ranges - CNS
!d4 = sum(ydot(first_surf_spec+n_surf_spec:first_surf_spec+n_surf_spec+n_surf_spec-1))
d4 = sum(ydot(first_bulk_spec:nspecies))


nml = aint(sum(y(first_surf_spec:nspecies))/nsites)
nml_s_array = aint(sum(s(first_surf_spec:nspecies)%abundance)/nsites)


IF (nml/=nml_old .AND. t/=told) THEN
    nml_int = aint(nml)
    PRINT *, "nml_old = ", nml_old
    PRINT *, "nml = ", nml
    PRINT *, "nml_s_array = ", nml_s_array
    PRINT *, "nml_int = ", nml_int
    PRINT *, "size of abundances_bulk = ",SIZE(abundances_bulk,1), " in dimension 1, and ",SIZE(abundances_bulk,2), " in dimension 2"
    PRINT *, "size of y = ",SIZE(y)
    IF (nml_int>0) abundances_bulk(1:nspecies,nml_int) = y(1:nspecies)
    IF (nml_int>0) timesteps_nml(nml_int) = t
    IF (nml>nml_max) nml_max = nml
!    CALL save_results_bulk
    write(57,'(4(1pe15.7))')t/3.155d7, nml, sum(y(first_surf_spec:nspecies))/nsites, nml_max
ENDIF
nml_old = nml

! If modeling an experiment and nml < 1, then end: not enough ice to study
IF ( (MODEL_EXPERIMENT.EQ.1).AND.(nml.LT.1) ) THEN
  PRINT *, "STOPPING SIMULATION! Number of Monolayers Less Than 1"
  CALL save_results_shingledecker
  CALL EXIT()
ENDIF

if (t/=told) then
  IF ( MODEL_EXPERIMENT .EQ. 0 ) THEN
    write(*,'(a11,1pe12.4,a5,1pe12.4,a5,1pe12.4)')'t (yrs.) = ',t/3.155e7 !,' d= ',gdens,' t = ',get_parameter(t/year, n_t_steps, time_t_array, t_array)
!    write(*,*) "Nml=",sum(y(first_surf_spec:nspecies))/nsites
!    write(*,*)'time, absolute coverage, alpha (fractional coverage), absolute coverage as sum(y(first_surf_spec:nspecies)), bulk abs. cov. by derivative, bult abs. cov. as sum of abundances'
!    write(*,'(7(1pe15.7))')t/3.155d7, dtran, diff_m2s_tot, sum(y(first_surf_spec:first_surf_spec+n_surf_spec-1)), dtran_fd, sum(y(first_surf_spec+n_surf_spec:first_surf_spec+n_surf_spec+n_surf_spec-1)), cov
!    write(*,*)"t/3.155d7, dtran, dtran_tmp, -diff_s2m_tot, diff_m2s_tot, dtran-dtran_tmp-diff_s2m_tot+diff_m2s_tot, d3, d4, dtran_chb"
!    write(*,'(10(1pe12.4,2x))')t/3.155d7, dtran, dtran_tmp, -diff_s2m_tot, diff_m2s_tot, dtran-dtran_tmp-diff_s2m_tot+diff_m2s_tot, d3, d4, dtran_chb
    PRINT *, "smol=",smol," and tot_surf_ab=",tot_surf_ab
    PRINT *, "bmol=",bmol," and tot_bulk_ab=",tot_bulk_ab
    IF ( ABS(smol-tot_surf_ab) .GT. EPSILON(smol) ) PRINT *, "!!!!!ERROR: bmol≠tot_bulk_ab, diff=",ABS(smol-tot_surf_ab)
    IF ( ABS(bmol-tot_bulk_ab) .GT. EPSILON(bmol) ) PRINT *, "!!!!!ERROR: bmol≠tot_bulk_ab, diff=",ABS(bmol-tot_bulk_ab)
    PRINT *, "dtran_freeze=",dtran_freeze
    PRINT *, "dtran_desorb=",dtran_desorb
    PRINT *, "dtran_chemdes=",dtran_chemdes
    PRINT *, "dtran=",dtran
    PRINT *, "diff_s2m_tot=",diff_s2m_tot
    PRINT *, "diff_m2s_tot=",diff_m2s_tot
!    write(*,*) "**********************************"
  ELSE
    1010 FORMAT(1X,A9,1X,ES12.4,A21,F12.4,A21)
!    write(*,'(a22,1pES12.4)') 'Fluence (ions/cm^2) = ',t*PHI_EXP
    write(*,*) "Nml=",sum(y(first_surf_spec:nspecies))/nsites
    write(*,'(a11,1pES12.4)') 'Time (s) = ',t
    IF ( MODEL_EXPERIMENT .EQ. 1 ) THEN
      PRINT '(A,ES12.4)', 'Fluence = ',t*PHI_EXP
    ENDIF
!    wrt = y(species_idx('bH2O      '))
!    wrt = (y(species_idx('gH2O      ')) + y(species_idx('bH2O      ')))
!    wrt = (s(species_idx('gO2       '))%abundance + s(species_idx('bO2       '))%abundance)
    wrt = initial_oxygen 
!     wrt = 1.0e20*1.0e-12
!    write(*,1010) "n(SO)   =",(y(species_idx('gSO       ')) + y(species_idx('bSO       ')))/1.0e20
!    write(*,1010) "n(SO2)  =",(y(species_idx('gSO2      ')) + y(species_idx('bSO2      ')))/1.0e20
!    write(*,1010) "n(SO3)  =",(y(species_idx('gSO3      ')) + y(species_idx('bSO3      ')))/1.0e20
!    write(*,1010) "n(S2O)  =",(y(species_idx('gS2O      ')) + y(species_idx('bS2O      ')))/1.0e20
!    WRITE(*,*) "-----GRAIN-----"
!    write(*,*)    "Current water = ", (100.0*wrt)/initial_water," % of initial"
!    write(*,1010) "n(H2O)grain  =",100*(y(species_idx('gH2O      ')) + y(species_idx('bH2O      ')))/wrt,"% wrt. total water : ",100.0*(y(species_idx('bH2O      '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(H2)grain   =",100*(y(species_idx('gH2       ')) + y(species_idx('bH2       ')))/wrt,"% wrt. total water : ",100.0*(y(species_idx('bH2       '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
    temp_atoms = 0
    temp_atoms_y = 0
    PRINT *, "In the normal species print-out location"
    DO i = first_surf_spec,nspecies
!      IF ( TRIM(s(i)%name) .NE. "e-" .AND. TRIM(s(i)%name) .NE. "ge-" .AND. TRIM(s(i)%name).NE."be-" ) THEN
        WRITE (*,'(A,A,ES9.2,A,F5.1,A,ES9.2,A,F5.1,A)') s(i)%name, "has an s-array abundance of ", s(i)%abundance, "(", ((s(i)%abundance*s(i)%natoms)/total_atoms)*100.0, &
                                                 "% init. O) and a y-array abundance of ",y(i), " (",((y(i)*s(i)%natoms)/total_atoms)*100.0,"% init. O)"
        temp_atoms = temp_atoms + s(i)%abundance*s(i)%natoms
        temp_atoms_y = temp_atoms_y + y(i)*s(i)%natoms
!      ENDIF
    ENDDO

    WRITE( *,'(A,F5.1,A)') "There are ", (temp_atoms/total_atoms)*100.0, "% of the initial atoms left in the s-array"
    WRITE( *,'(A,F8.4,A)') "There are ", (temp_atoms_y/total_atoms)*100.0000, "% of the initial atoms left in the y-array"

    OPEN(53066, file='y-array_over_time.csv', status='old', position='append') ! , iostat = err_status, iomsg=err_iomsg

    ! IF(err_status /= 0) THEN
    !   WRITE (*,*) 'Error ', TRIM(err_iomsg)
    !   STOP
    ! END IF

    WRITE(53066, '(1ES12.4,A,F8.4)') t, ", ", (temp_atoms_y/total_atoms)*100.0000
    CLOSE(53066)    

!    write(*,1010) "n(O2)grain     =",100*(y(species_idx('gO2       ')) + y(species_idx('bO2       ')))/wrt,"% wrt. initial O2 : " ,100.0*(y(species_idx('bO2       '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(O2*)grain    =",100*(y(species_idx('gO2*      ')) + y(species_idx('bO2*      ')))/wrt,"% wrt. initial O2 : " ,100.0*(y(species_idx('bO2*      '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(O2+)grain    =",100*(y(species_idx('gO2+      ')) + y(species_idx('bO2+      ')))/wrt,"% wrt. initial O2 : " ,100.0*(y(species_idx('bO2+      '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(O2-)grain    =",100*(y(species_idx('gO2-      ')) + y(species_idx('bO2-       ')))/wrt,"% wrt. initial O2 : ",100.0*(y(species_idx('bO2-      '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(O2)grain   =",100*(y(species_idx('gO2       ')) + y(species_idx('bO2       ')))/wrt,"% wrt. total water : ",100.0*(y(species_idx('bO2       '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(O3)grain     =",100*(y(species_idx('gO3       ')) +  y(species_idx('bO3       ')))/wrt,"% wrt. initial O2 : ",100.0*(y(species_idx('bO3       '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(O3*)grain    =",100*(y(species_idx('gO3*      ')) +  y(species_idx('bO3*      ')))/wrt,"% wrt. initial O2 : ",100.0*(y(species_idx('bO3*      '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(O3+)grain    =",100*(y(species_idx('gO3+      ')) + y(species_idx('bO3+      ')))/wrt,"% wrt. initial O2 : ",100.0*(y(species_idx('bO3+      '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(O3-)grain    =",100*(y(species_idx('gO3-      ')) + y(species_idx('bO3-      ')))/wrt,"% wrt. initial O2 : ",100.0*(y(species_idx('bO3-      '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(O3)grain   =",100*(y(species_idx('gO3       ')) + y(species_idx('bO3       ')))/wrt,"% wrt. total water : ",100.0*(y(species_idx('bO3       '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(H2O2)grain =",100*(y(species_idx('gH2O2     ')) + y(species_idx('bH2O2     ')))/wrt,"% wrt. total water : ",100.0*(y(species_idx('bH2O2     '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(H)grain    =",100*(y(species_idx('gH        ')) + y(species_idx('bH        ')))/wrt,"% wrt. total water : ",100.0*(y(species_idx('bH        '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(O)grain      =",100*(y(species_idx('gO        ')) + y(species_idx('bO        ')))/wrt,"% wrt. initial O2 : ",100.0*(y(species_idx('bO        '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(O*)grain     =",100*(y(species_idx('gO*       ')) + y(species_idx('bO*       ')))/wrt,"% wrt. initial O2 : ",100.0*(y(species_idx('bO*       '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(O+)grain     =",100*(y(species_idx('gO+       ')) + y(species_idx('bO+       ')))/wrt,"% wrt. initial O2 : ",100.0*(y(species_idx('bO+       '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(O-)grain     =",100*(y(species_idx('gO-       ')) + y(species_idx('bO-       ')))/wrt,"% wrt. initial O2 : ",100.0*(y(species_idx('bO-       '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(e-)grain     =",100*(y(species_idx('ge-       ')) + y(species_idx('be-       ')))/wrt,"% wrt. initial O2 : ",100.0*(y(species_idx('be-       '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(O)grain    =",100*(y(species_idx('gO        ')) + y(species_idx('bO        ')))/wrt,"% wrt. total water : ",100.0*(y(species_idx('bO        '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(OH)grain   =",100*(y(species_idx('gOH       ')) + y(species_idx('bOH       ')))/wrt,"% wrt. total water : ",100.0*(y(species_idx('bOH       '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(HO2)grain  =",100*(y(species_idx('gHO2      ')) + y(species_idx('bHO2      ')))/wrt,"% wrt. total water : ",100.0*(y(species_idx('bHO2      '))/(tot_surf_ab + tot_bulk_ab)),"% wrt. total ice ab."
!    write(*,1010) "n(HO3)grain  =",100*(y(species_idx('gHO3      ')) + y(species_idx('bHO3      ')))/wrt
!    WRITE(*,*) "-----GAS-----"
!    write(*,1010) "n(H)gas      =",100*y(species_idx('H         '))/wrt
!    write(*,1010) "n(H2)gas     =",100*y(species_idx('H2        '))/wrt
!    write(*,1010) "n(O)gas      =",100*y(species_idx('O         '))/wrt
!    write(*,1010) "n(O+)gas      =",100*y(species_idx('O+        '))/wrt
!    write(*,1010) "n(O-)gas      =",100*y(species_idx('O-        '))/wrt
!    write(*,1010) "n(O2)gas     =",100*y(species_idx('O2        '))/wrt
!    write(*,1010) "n(O2+)gas     =",100*y(species_idx('O2+       '))/wrt
!    write(*,1010) "n(O2-)gas     =",100*y(species_idx('O2-       '))/wrt
!    write(*,1010) "n(O3)gas     =",100*y(species_idx('O3        '))/wrt
!    write(*,1010) "n(O3+)gas     =",100*y(species_idx('O3+       '))/wrt
!    write(*,1010) "n(O3-)gas     =",100*y(species_idx('O3-       '))/wrt
!    write(*,1010) "n(OH)gas     =",100*y(species_idx('OH        '))/wrt
!    write(*,1010) "n(H2O)gas    =",100*y(species_idx('H2O       '))/wrt
!    write(*,1010) "n(HO2)gas    =",100*y(species_idx('HO2       '))/wrt
!    write(*,1010) "n(HO3)gas    =", y(species_idx('HO3       '))/wrt
!    write(*,1010) "n(HS)   =",(y(species_idx('gHS       ')) + y(species_idx('bHS       ')))/1.0e20
!    write(*,1010) "n(H2S)  =",(y(species_idx('gH2S      ')) + y(species_idx('bH2S      ')))/1.0e20
!    write(*,1010) "n(S2H)  =",(y(species_idx('gS2H      ')) + y(species_idx('bS2H      ')))/1.0e20
!    write(*,1010) "n(H2S2) =",(y(species_idx('gH2S2     ')) + y(species_idx('bH2S2     ')))/1.0e20
!    write(*,*) "---------------------------------"
!    write(*,1010) "n(S)    =",(y(species_idx('gS        ')) + y(species_idx('bS        ')))/1.0e20
!    write(*,1010) "n(S2)   =",(y(species_idx('gS2       ')) + y(species_idx('bS2       ')))/1.0e20
!    write(*,1010) "n(S3)   =",(y(species_idx('gS3       ')) + y(species_idx('bS3       ')))/1.0e20
!    write(*,1010) "n(S4)   =",(y(species_idx('gS4       ')) + y(species_idx('bS4       ')))/1.0e20
!    write(*,1010) "n(S5)   =",(y(species_idx('gS5       ')) + y(species_idx('bS5       ')))/1.0e20
!    write(*,1010) "n(S6)   =",(y(species_idx('gS6       ')) + y(species_idx('bS6       ')))/1.0e20
!    write(*,1010) "n(S7)   =",(y(species_idx('gS7       ')) + y(species_idx('bS7       ')))/1.0e20
!    write(*,1010) "n(S8)   =",(y(species_idx('gS8       ')) + y(species_idx('bS8       ')))/1.0e20
    write(*,*) "**********************************"
    write(*,*) "**********************************"
  ENDIF
  write(41,*)t/3.155d7, dtran
  write(42,'(7(1pe15.7))')t/3.155d7, dtran, diff_m2s_tot, sum(y(first_surf_spec:first_surf_spec+n_surf_spec-1)), dtran_fd, sum(y(first_surf_spec+n_surf_spec:first_surf_spec+n_surf_spec+n_surf_spec-1)), cov
  write(43,*)t/3.155d7, transit
  write(44,'(10(1pe12.4,2x))')t/3.155d7, dtran, dtran_tmp, -diff_s2m_tot, diff_m2s_tot, dtran-dtran_tmp-diff_s2m_tot+diff_m2s_tot, d3, d4, dtran_chb
  IF ( MODEL_EXPERIMENT .EQ. 0 ) THEN
    write(55,'(3(1pe12.4,2x))')t/3.155d7, sum(y(first_surf_spec:nspecies))/nsites, y(species_idx('gH2O      '))/sum(y(first_surf_spec:first_surf_spec+n_surf_spec-1))
  ENDIF
  told=t
endif

RETURN
END SUBROUTINE re

SUBROUTINE sparse_jac(N, T, Y, IA, JA, NZ, P)
IMPLICIT NONE
INTEGER                               :: N, NZ, i, j, k, l
INTEGER, DIMENSION(*)                 :: IA, JA
REAL(wp)                              :: T
REAL(wp), DIMENSION(N)                :: Y, ddtran_dxj, dtransit_dxj
REAL(wp), DIMENSION(*)                :: P
REAL(wp), DIMENSION(:,:), ALLOCATABLE :: pd

ALLOCATE(pd(nspecies,nspecies))
pd(:,:) = 0.0d0

DO j = 1, nreactions
  IF (r(j)%ir2/=0) THEN
    PD(r(j)%ir1,r(j)%ir1) = PD(r(j)%ir1,r(j)%ir1) - r(j)%rate*y(r(j)%ir2)
    PD(r(j)%ir2,r(j)%ir1) = PD(r(j)%ir2,r(j)%ir1) - r(j)%rate*y(r(j)%ir2)
    PD(r(j)%ip1,r(j)%ir1) = PD(r(j)%ip1,r(j)%ir1) + r(j)%rate*y(r(j)%ir2)
    IF (r(j)%ip2/=0) PD(r(j)%ip2,r(j)%ir1) = PD(r(j)%ip2,r(j)%ir1) + r(j)%rate*y(r(j)%ir2)
    IF (r(j)%ip3/=0) PD(r(j)%ip3,r(j)%ir1) = PD(r(j)%ip3,r(j)%ir1) + r(j)%rate*y(r(j)%ir2)
    IF (r(j)%ip4/=0) PD(r(j)%ip4,r(j)%ir1) = PD(r(j)%ip4,r(j)%ir1) + r(j)%rate*y(r(j)%ir2)
    IF (r(j)%ip5/=0) PD(r(j)%ip5,r(j)%ir1) = PD(r(j)%ip5,r(j)%ir1) + r(j)%rate*y(r(j)%ir2)

    PD(r(j)%ir1,r(j)%ir2) = PD(r(j)%ir1,r(j)%ir2) - r(j)%rate*y(r(j)%ir1)
    PD(r(j)%ir2,r(j)%ir2) = PD(r(j)%ir2,r(j)%ir2) - r(j)%rate*y(r(j)%ir1)
    PD(r(j)%ip1,r(j)%ir2) = PD(r(j)%ip1,r(j)%ir2) + r(j)%rate*y(r(j)%ir1)
    IF (r(j)%ip2/=0) PD(r(j)%ip2,r(j)%ir2) = PD(r(j)%ip2,r(j)%ir2) + r(j)%rate*y(r(j)%ir1)
    IF (r(j)%ip3/=0) PD(r(j)%ip3,r(j)%ir2) = PD(r(j)%ip3,r(j)%ir2) + r(j)%rate*y(r(j)%ir1)
    IF (r(j)%ip4/=0) PD(r(j)%ip4,r(j)%ir2) = PD(r(j)%ip4,r(j)%ir2) + r(j)%rate*y(r(j)%ir1)
    IF (r(j)%ip5/=0) PD(r(j)%ip5,r(j)%ir2) = PD(r(j)%ip5,r(j)%ir2) + r(j)%rate*y(r(j)%ir1)
  ELSE
    PD(r(j)%ir1,r(j)%ir1) = PD(r(j)%ir1,r(j)%ir1) - r(j)%rate
    PD(r(j)%ip1,r(j)%ir1) = PD(r(j)%ip1,r(j)%ir1) + r(j)%rate
    IF (r(j)%ip2/=0) PD(r(j)%ip2,r(j)%ir1) = PD(r(j)%ip2,r(j)%ir1) + r(j)%rate
    IF (r(j)%ip3/=0) PD(r(j)%ip3,r(j)%ir1) = PD(r(j)%ip3,r(j)%ir1) + r(j)%rate
    IF (r(j)%ip4/=0) PD(r(j)%ip4,r(j)%ir1) = PD(r(j)%ip4,r(j)%ir1) + r(j)%rate
    IF (r(j)%ip5/=0) PD(r(j)%ip5,r(j)%ir1) = PD(r(j)%ip5,r(j)%ir1) + r(j)%rate
  ENDIF
ENDDO

ddtran_dxj(1:N) = 0.0d0
dtransit_dxj(1:N) = 0.0d0

DO j = 1, N
    ddtran_dxj(j) = sum(PD(first_surf_spec:first_surf_spec+n_surf_spec-1,j))
ENDDO

if (tdust>150.0d0) then
    dtran = 0.0d0
    ddtran_dxj(1:n) = 0.0d0
endif

IF (dtran>0.0d0) THEN

    DO i = first_surf_spec, first_surf_spec + n_surf_spec - 1
        DO j = 1, nspecies
            PD(i,j) = PD (i,j) - (ddtran_dxj(j)*y(i)+dtran*delta(i,j))/(nsites*n_s_ml)
            PD(i+n_surf_spec,j) = PD(i+n_surf_spec,j) + (ddtran_dxj(j)*y(i)+dtran*delta(i,j))/(nsites*n_s_ml)
        ENDDO
    ENDDO

ENDIF

IF (dtran<0.0d0) THEN

    DO i = first_surf_spec, first_surf_spec + n_surf_spec - 1
        DO j = 1, nspecies

            IF (bmol<smol) THEN
                IF (j>=first_surf_spec .AND. j<=first_surf_spec+n_surf_spec-1) THEN
                    PD(i,j) = PD(i,j) - ddtran_dxj(j)*y(i+n_surf_spec)/smol+dtran*(delta(i+n_surf_spec,j)/smol-y(i+n_surf_spec)/smol**2.0d0)
                    PD(i+n_surf_spec,j) = PD(i+n_surf_spec,j) + ddtran_dxj(j)*y(i+n_surf_spec)/smol+dtran*(delta(i+n_surf_spec,j)/smol-y(i+n_surf_spec)/smol**2.0d0)
                ELSE
                    PD(i,j) = PD(i,j) - ddtran_dxj(j)*y(i+n_surf_spec)/smol+dtran*(delta(i+n_surf_spec,j)/smol-0.0d0)
                    PD(i+n_surf_spec,j) = PD(i+n_surf_spec,j) + ddtran_dxj(j)*y(i+n_surf_spec)/smol+dtran*(delta(i+n_surf_spec,j)/smol-0.0d0)
                ENDIF
            ELSE
                IF (j>=first_surf_spec+n_surf_spec) THEN
                    PD(i,j) = PD(i,j) - ddtran_dxj(j)*y(i+n_surf_spec)/bmol+dtran*(delta(i+n_surf_spec,j)/bmol-y(i+n_surf_spec)/(bmol**2.0d0))
                    PD(i+n_surf_spec,j) = PD(i+n_surf_spec,j) + ddtran_dxj(j)*y(i+n_surf_spec)/bmol+dtran*(delta(i+n_surf_spec,j)/bmol-y(i+n_surf_spec)/(bmol**2.0d0))
                ELSE
                    PD(i,j) = PD(i,j) - ddtran_dxj(j)*y(i+n_surf_spec)/bmol+dtran*delta(i+n_surf_spec,j)/bmol
                    PD(i+n_surf_spec,j) = PD(i+n_surf_spec,j) + ddtran_dxj(j)*y(i+n_surf_spec)/bmol+dtran*delta(i+n_surf_spec,j)/bmol
                ENDIF
            ENDIF

        ENDDO
    ENDDO

ENDIF

k = 0

DO i = 1, n
  DO j = 1, n
    IF (i==j .OR. PD(i,j)/=0.0d0) k = k + 1
  ENDDO
ENDDO
!print*, k, nz
IF (NZ==0) THEN
  NZ = 250000 !k
ELSE
  IA(1) = 1
  k = 0
  l = 0

  DO i = 1, n
    l = 0
    DO j = 1, n
      IF (i==j .OR. PD(j,i)/=0.0d0) THEN
        k = k + 1
        l = l + 1
        JA(k) = j
        P(k) = PD(j,i)
      ENDIF
    ENDDO
    IA(i+1) = IA(i) + l
  ENDDO
ENDIF

DEALLOCATE(pd)

END SUBROUTINE sparse_jac

SUBROUTINE mre(neq, t, y, ydot)
INTEGER                  :: neq, i, j
REAL(wp)                 :: t, rr, told, Rmod, Racc1, Racc2, fAB, cAB, dAB1, dAB2, eAB1, eAB2, fAB1, fAB2, eta1, eta2
REAL(wp), DIMENSION(neq) :: y, ydot

ydot(:) = 0.0d0
IF (delta_rho==1 .OR. delta_t==1) CALL calc_rates(t)

DO i = 1,nreactions
  rr = 0.0d0
  IF (r(i)%ir2==0) THEN
    rr = r(i)%rate*y(r(i)%ir1)
  ELSE
    rr = r(i)%rate*y(r(i)%ir1)*y(r(i)%ir2)
  ENDIF

    IF (r(i)%rtype==13) THEN

    Racc1=s(s(r(i)%ir1)%gas_idx)%racc*y(s(r(i)%ir1)%gas_idx)*r(i)%rate/(r(i)%rate+s(r(i)%ir1)%rdes+s(r(i)%ir2)%rdes)
    Racc2=s(s(r(i)%ir2)%gas_idx)%racc*y(s(r(i)%ir2)%gas_idx)*r(i)%rate/(r(i)%rate+s(r(i)%ir1)%rdes+s(r(i)%ir2)%rdes)

    DO j=first_surfreact, nreactions

    IF (r(i)%ir1==r(j)%ip1 .OR. r(i)%ir1==r(j)%ip2 .OR. r(i)%ir1==r(j)%ip3) THEN

      fAB=100.0d0/(100.0d0+(y(r(j)%ir1)+y(r(j)%ir2))**6)
      Racc1=Racc1+(1.0d0-fAB)*r(j)%rate*y(r(j)%ir1)*y(r(j)%ir2)

    ENDIF

    IF (r(i)%ir2==r(j)%ip1 .OR. r(i)%ir2==r(j)%ip2 .OR. r(i)%ir2==r(j)%ip3) THEN

      fAB=100.0d0/(100.0d0+(y(r(j)%ir1)+y(r(j)%ir2))**6)
      Racc2=Racc2+(1.0d0-fAB)*r(j)%rate*y(r(j)%ir1)*y(r(j)%ir2)

    ENDIF

    ENDDO

    cAB=r(i)%rate/(r(i)%rate+s(r(i)%ir1)%rdes+s(r(i)%ir2)%rdes)

    dAB1=s(r(i)%ir2)%rdes/(r(i)%rate+s(r(i)%ir1)%rdes+s(r(i)%ir2)%rdes)
    dAB2=s(r(i)%ir1)%rdes/(r(i)%rate+s(r(i)%ir1)%rdes+s(r(i)%ir2)%rdes)

    eAB1=Racc2/(Racc2+s(r(i)%ir1)%rdes)
    eAB2=Racc1/(Racc1+s(r(i)%ir2)%rdes)

    fAB1=dAB1*eAB1
    fAB2=dAB2*eAB2

    eta1=1.0d0 !cAB !+fAB1/2.0d0+fAB**2/3.0d0+fAB**3/4.0d0
    eta2=1.0d0 !cAB !+fAB1/2.0d0+fAB**2/3.0d0+fAB**3/4.0d0

    IF (r(i)%ir1.ne.r(i)%ir2) THEN
      Rmod=Racc2*y(r(i)%ir1)*eta1+Racc1*y(r(i)%ir2)*eta2
    ELSE
      Rmod=Racc1*y(r(i)%ir1)*eta1
    ENDIF

    fAB=100.0d0/(100.0d0+(y(r(i)%ir1)+y(r(i)%ir2))**6)

    IF (Rmod/=0.0d0 .AND. rr/=0.0d0) &
           rr=(1.0d0/(1.0d0+(rr/Rmod)**2)*rr+ &
             (1.0d0-1.0d0/(1.0d0+(rr/Rmod)**2))* &
               Rmod*fAB+r(i)%rate*y(r(i)%ir1)*y(r(i)%ir2)*(1.0d0-fAB))*r(i)%alpha

    ENDIF

  mre_terms(i) = rr

  ydot(r(i)%ir1) = ydot(r(i)%ir1) - rr
  IF (r(i)%ir2/=0) ydot(r(i)%ir2) = ydot(r(i)%ir2) - rr
  ydot(r(i)%ip1) = ydot(r(i)%ip1) + rr
  IF (r(i)%ip2/=0) ydot(r(i)%ip2) = ydot(r(i)%ip2) + rr
  IF (r(i)%ip3/=0) ydot(r(i)%ip3) = ydot(r(i)%ip3) + rr
  IF (r(i)%ip4/=0) ydot(r(i)%ip4) = ydot(r(i)%ip4) + rr
  IF (r(i)%ip5/=0) ydot(r(i)%ip5) = ydot(r(i)%ip5) + rr
ENDDO

if (t/=told) then
  print*, 't (yrs.) = ',t/3.155e7
  told=t
endif
END SUBROUTINE mre

SUBROUTINE open_analytics_files
IMPLICIT NONE
INTEGER :: i

DO i = 1, n_det_spec
  OPEN(i+77, FILE='analytics_'//s_det_study(i)%name(1:LEN_TRIM(s_det_study(i)%name)), STATUS='UNKNOWN', ACCESS='APPEND')
  REWIND(i+77)
ENDDO

END SUBROUTINE open_analytics_files

SUBROUTINE close_analytics_files
IMPLICIT NONE
INTEGER :: i

DO i = 1, n_det_spec
  CLOSE(i+77)
ENDDO

END SUBROUTINE close_analytics_files

SUBROUTINE save_analytics(neq,y,t, nt)
IMPLICIT NONE
INTEGER                  :: i, neq, nt
REAL(wp), DIMENSION(neq) :: y
REAL(wp)                 :: t

DO i = 1, n_det_spec
  CALL get_analytics_on_species(s_det_study(i)%name,neq,y,t,i+77, nt, i)
ENDDO

END SUBROUTINE save_analytics

SUBROUTINE get_analytics_on_species(s, neq, y, t, fileid, ntt, ndss)
IMPLICIT NONE
CHARACTER*10 s

REAL(wp)                   :: rform, rdest, rtmp, rr, t, b_r_divider
REAL(wp), DIMENSION(10000) :: rates, rates_s
INTEGER, DIMENSION(10000)  :: irates, irates_s
INTEGER                    :: j,c, l, m, itmp, itmp2, neq, fileid, ntt, ndss
REAL(wp), DIMENSION(neq)   :: y

c = 1
rform = 0.0d0
rdest = 0.0d0

DO j = 1, nreactions
  rr = 0.0d0
  IF (r(j)%ir2==0 .OR. r(j)%ir2==17 .OR. r(j)%ir2==18) THEN
    rr = r(j)%rate*y(r(j)%ir1)
  ELSE
    rr = r(j)%rate*y(r(j)%ir1)*y(r(j)%ir2)
  ENDIF

    !Bulk reactions:
    IF (r(j)%rtype==14 .OR. r(j)%rtype==16) THEN
        IF (bmol<nsites) THEN
                rr = rr/nsites
        ELSE
                rr = rr/bmol
        ENDIF
    ENDIF

  IF (eqtype==2) rr = mre_terms(j)

  IF (des_reactive_type==2 .AND. (r(j)%rtype==13 .OR. r(j)%rtype==15)) THEN
      IF (r(j)%p1==s .AND. r(j)%p1(1:1)=='g') rr = rr - rd_v2_terms(1,j)
      IF (r(j)%p2==s .AND. r(j)%p2(1:1)=='g') rr = rr - rd_v2_terms(2,j)
      IF (r(j)%p3==s .AND. r(j)%p3(1:1)=='g') rr = rr - rd_v2_terms(3,j)
      IF (r(j)%p4==s .AND. r(j)%p4(1:1)=='g') rr = rr - rd_v2_terms(4,j)
      IF (r(j)%p5==s .AND. r(j)%p5(1:1)=='g') rr = rr - rd_v2_terms(5,j)

      IF (r(j)%p1==s .AND. r(j)%p1(1:1)/='g') rr = rr + rd_v2_terms(1,j)
      IF (r(j)%p2==s .AND. r(j)%p2(1:1)/='g') rr = rr + rd_v2_terms(2,j)
      IF (r(j)%p3==s .AND. r(j)%p3(1:1)/='g') rr = rr + rd_v2_terms(3,j)
      IF (r(j)%p4==s .AND. r(j)%p4(1:1)/='g') rr = rr + rd_v2_terms(4,j)
      IF (r(j)%p5==s .AND. r(j)%p5(1:1)/='g') rr = rr + rd_v2_terms(5,j)
  ENDIF

  ! Modified by C. N. Shingledecker
!  IF ( ( r(j)%rtype .NE. 15 ) .AND. ( r(j)%rtype .NE. 16 ) .AND. ( r(j)%rtype .NE. 17 ) .AND. ( r(j)%rtype .NE. 18 ) )THEN
    IF (r(j)%r1==s .OR. r(j)%r2==s .OR. r(j)%p1==s .OR. r(j)%p2==s .OR. r(j)%p3==s .OR. r(j)%p4==s .OR. r(j)%p5==s) THEN
      IF (r(j)%r1==s .OR. r(j)%r2==s) THEN
        rdest = rdest + rr
      ELSE
        rform = rform + rr
      ENDIF
      rates(c) = rr
      irates(c) = j
      c = c + 1
    ENDIF
!  ENDIF
ENDDO

DO l = 1, c
  rtmp = 0.0d0
  itmp = 0
  DO m = 1, c
    IF (rates(m)>rtmp) THEN
      rtmp = rates(m)
      itmp = irates(m)
      itmp2 = m
    ENDIF
  ENDDO
  IF (itmp/=0) THEN
    rates_s(l) = rtmp
    irates_s(l) = itmp

    reaction_importance(ndss, irates_s(l), ntt) = rates_s(l)/(rdest+rform)

        !Bulk reactions:
        b_r_divider = 1.0d0
        IF (r(irates_s(l))%rtype==14 .OR. r(irates_s(l))%rtype==16) THEN
            IF (bmol<nsites) THEN
                b_r_divider = nsites
            ELSE
                b_r_divider = bmol
            ENDIF
        ENDIF

    rates(itmp2) = 0.0d0


    ! NB: This overwrite will make the bulk rates display incorrectly
    IF ( r(irates_s(l))%rtype .EQ. 1 ) THEN
      b_r_divider = ddens
    ENDIF

    IF (r(irates_s(l))%r1==s .OR. r(irates_s(l))%r2==s) THEN
      WRITE(fileid,'(i4,1x,i4,1x, 7a10, 1x, es10.4, 1x, f6.2, 3x, i2, 3x, es10.4, 3x, es10.4)')l, r(irates_s(l))%idx, r(irates_s(l))%r1, r(irates_s(l))%r2, r(irates_s(l))%p1, r(irates_s(l))%p2, r(irates_s(l))%p3, r(irates_s(l))%p4, r(irates_s(l))%p5, rates_s(l), -rates_s(l)/(rdest+rform)*100., r(irates_s(l))%rtype, r(irates_s(l))%rate/b_r_divider, b_r_divider
!      IF (ABS(rates_s(l)/(rdest+rform)*100.)>5.0d-3) WRITE(fileid,'(i4,1x,i4,1x, 7a10, 1x, e10.4, 1x, f6.2, 3x, i2, 3x, e10.4, 3x, e10.4)')l, r(irates_s(l))%idx, r(irates_s(l))%r1, r(irates_s(l))%r2, r(irates_s(l))%p1, r(irates_s(l))%p2, r(irates_s(l))%p3, r(irates_s(l))%p4, r(irates_s(l))%p5, rates_s(l), -rates_s(l)/(rdest+rform)*100., r(irates_s(l))%rtype, r(irates_s(l))%rate/b_r_divider, b_r_divider
    ELSE
      WRITE(fileid,'(i4,1x,i4,1x, 7a10, 1x, es10.4, 1x, f6.2, 3x, i2, 3x, es10.4, 3x, es10.4)')l, r(irates_s(l))%idx, r(irates_s(l))%r1, r(irates_s(l))%r2, r(irates_s(l))%p1, r(irates_s(l))%p2, r(irates_s(l))%p3, r(irates_s(l))%p4, r(irates_s(l))%p5, rates_s(l), rates_s(l)/(rdest+rform)*100., r(irates_s(l))%rtype, r(irates_s(l))%rate/b_r_divider, b_r_divider
!      IF (ABS(rates_s(l)/(rdest+rform)*100.)>5.0d-3) WRITE(fileid,'(i4,1x,i4,1x, 7a10, 1x, e10.4, 1x, f6.2, 3x, i2, 3x, e10.4, 3x, e10.4)')l, r(irates_s(l))%idx, r(irates_s(l))%r1, r(irates_s(l))%r2, r(irates_s(l))%p1, r(irates_s(l))%p2, r(irates_s(l))%p3, r(irates_s(l))%p4, r(irates_s(l))%p5, rates_s(l), rates_s(l)/(rdest+rform)*100., r(irates_s(l))%rtype, r(irates_s(l))%rate/b_r_divider, b_r_divider
    ENDIF
  ENDIF
ENDDO

IF ( MODEL_EXPERIMENT .EQ. 1 ) THEN
  WRITE(fileid,'(a24,E10.4)') "============ FLUENCE =  ",(t*PHI_EXP)/1.0e14
ELSE
  WRITE(fileid,'(a24,e10.4)') "============ TIME (yr) =",t/year
ENDIF


END SUBROUTINE get_analytics_on_species

SUBROUTINE jac(NEQ,T,Y,ML,MU,PD,NRPD)
IMPLICIT NONE
INTEGER NEQ, ML, MU, NRPD, j
REAL(wp) PD, T, Y
DIMENSION Y(NEQ), PD(NRPD,NEQ)

PD(:,:) = 0.0d0

DO j = 1, nreactions
  IF (r(j)%ir2/=0) THEN
    PD(r(j)%ir1,r(j)%ir1) = PD(r(j)%ir1,r(j)%ir1) - r(j)%rate*y(r(j)%ir2)
    PD(r(j)%ir2,r(j)%ir1) = PD(r(j)%ir2,r(j)%ir1) - r(j)%rate*y(r(j)%ir2)
    PD(r(j)%ip1,r(j)%ir1) = PD(r(j)%ip1,r(j)%ir1) + r(j)%rate*y(r(j)%ir2)
    IF (r(j)%ip2/=0) PD(r(j)%ip2,r(j)%ir1) = PD(r(j)%ip2,r(j)%ir1) + r(j)%rate*y(r(j)%ir2)
    IF (r(j)%ip3/=0) PD(r(j)%ip3,r(j)%ir1) = PD(r(j)%ip3,r(j)%ir1) + r(j)%rate*y(r(j)%ir2)
    IF (r(j)%ip4/=0) PD(r(j)%ip4,r(j)%ir1) = PD(r(j)%ip4,r(j)%ir1) + r(j)%rate*y(r(j)%ir2)
    IF (r(j)%ip5/=0) PD(r(j)%ip5,r(j)%ir1) = PD(r(j)%ip5,r(j)%ir1) + r(j)%rate*y(r(j)%ir2)

    PD(r(j)%ir1,r(j)%ir2) = PD(r(j)%ir1,r(j)%ir2) - r(j)%rate*y(r(j)%ir1)
    PD(r(j)%ir2,r(j)%ir2) = PD(r(j)%ir2,r(j)%ir2) - r(j)%rate*y(r(j)%ir1)
    PD(r(j)%ip1,r(j)%ir2) = PD(r(j)%ip1,r(j)%ir2) + r(j)%rate*y(r(j)%ir1)
    IF (r(j)%ip2/=0) PD(r(j)%ip2,r(j)%ir2) = PD(r(j)%ip2,r(j)%ir2) + r(j)%rate*y(r(j)%ir1)
    IF (r(j)%ip3/=0) PD(r(j)%ip3,r(j)%ir2) = PD(r(j)%ip3,r(j)%ir2) + r(j)%rate*y(r(j)%ir1)
    IF (r(j)%ip4/=0) PD(r(j)%ip4,r(j)%ir2) = PD(r(j)%ip4,r(j)%ir2) + r(j)%rate*y(r(j)%ir1)
    IF (r(j)%ip5/=0) PD(r(j)%ip5,r(j)%ir2) = PD(r(j)%ip5,r(j)%ir2) + r(j)%rate*y(r(j)%ir1)
  ELSE
    PD(r(j)%ir1,r(j)%ir1) = PD(r(j)%ir1,r(j)%ir1) - r(j)%rate
    PD(r(j)%ip1,r(j)%ir1) = PD(r(j)%ip1,r(j)%ir1) + r(j)%rate
    IF (r(j)%ip2/=0) PD(r(j)%ip2,r(j)%ir1) = PD(r(j)%ip2,r(j)%ir1) + r(j)%rate
    IF (r(j)%ip3/=0) PD(r(j)%ip3,r(j)%ir1) = PD(r(j)%ip3,r(j)%ir1) + r(j)%rate
    IF (r(j)%ip4/=0) PD(r(j)%ip4,r(j)%ir1) = PD(r(j)%ip4,r(j)%ir1) + r(j)%rate
    IF (r(j)%ip5/=0) PD(r(j)%ip5,r(j)%ir1) = PD(r(j)%ip5,r(j)%ir1) + r(j)%rate
  ENDIF
ENDDO

RETURN
END SUBROUTINE jac

END
