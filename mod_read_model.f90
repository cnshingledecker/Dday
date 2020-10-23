MODULE read_model
USE global_variables
USE read_rate06
!The module contains the subroutine which reads the model.inp file with the model set up
IMPLICIT NONE

CONTAINS

SUBROUTINE read_model_setup
IMPLICIT NONE
INTEGER      :: i, i_e
INTEGER      :: Ngas,Nsurf,Nbulk
INTEGER      :: io
CHARACTER*23 :: keyname
CHARACTER*6  :: skeyvalue
CHARACTER*40 :: keycomment
CHARACTER*13 :: tmp
REAL(wp)     :: positive_ion_abundance, negative_ion_and_electron_abundance, electron_correction

OPEN(1,FILE='model.inp',STATUS='OLD',ERR=100)
READ(1,*)
READ(1,*)
READ(1,1000)keyname, sing_mult, keycomment
READ(1,1000)keyname, is_disk_model, keycomment
READ(1,1000)keyname, RADIOLYSIS, keycomment
READ(1,1000)keyname, SUPRATHERMAL, keycomment
READ(1,1000)keyname, delta_rho, keycomment
READ(1,1000)keyname, delta_t, keycomment
READ(1,1000)keyname, delta_tdust, keycomment
READ(1,1000)keyname, delta_g0, keycomment
READ(1,1000)keyname, delta_avst, keycomment
READ(1,1000)keyname, delta_avis, keycomment
READ(1,1000)keyname, delta_zetacr, keycomment
READ(1,1000)keyname, delta_zetax, keycomment
READ(1,1000)keyname, delta_selfshield, keycomment
READ(1,1000)keyname, MODEL_EXPERIMENT, keycomment
READ(1,1000)keyname, FAST_BULK, keycomment
READ(1,1000)keyname, FAST_ATOMS, keycomment
READ(1,1000)keyname, PHOTOEXC, keycomment
READ(1,1000)keyname, PHOTOION, keycomment
READ(1,1000)keyname, FIXED_DVAL, keycomment
READ(1,1000)keyname, FIXED_NU, keycomment
READ(1,1000)keyname, DISABLE_DESORB, keycomment

READ(1,*)
READ(1,1002)keyname, chem_file
READ(1,1001)keyname, rtol, keycomment
READ(1,1001)keyname, atol, keycomment
READ(1,*)
READ(1,1001)keyname, n_s_ml
READ(1,*)
READ(1,1001)keyname, DVAL, keycomment
READ(1,1001)keyname, ICE_THICK, keycomment
READ(1,1001)keyname, RHO_ICE, keycomment
READ(1,1001)keyname, gdens, keycomment
READ(1,1001)keyname, t, keycomment
READ(1,1001)keyname, tdust, keycomment
READ(1,1001)keyname, G0_stellar, keycomment
IF (is_disk_model == 0) G0_stellar = 0.0d0
READ(1,1001)keyname, AvSt, keycomment
IF (is_disk_model == 0) AvSt = 0.0d0
READ(1,1001)keyname, AvIS, keycomment
READ(1,1001)keyname, ZetaCR, keycomment
ZetaCR = ZetaCR + 6.5D-19
READ(1,1001)keyname, ZetaX, keycomment
IF (is_disk_model == 0) ZetaX = 0.0d0
READ(1,1001)keyname, PHI_EXP, keycomment
READ(1,1001)keyname, SE_EXP, keycomment
READ(1,1001)keyname, TRIAL_NU, keycomment
READ(1,1000)keyname, des_t, keycomment
READ(1,1000)keyname, des_crp, keycomment
READ(1,1000)keyname, des_photon, keycomment
READ(1,1001)keyname, ph_yield, keycomment
READ(1,1001)keyname, des_reactive, keycomment
READ(1,1000)keyname, des_reactive_type, keycomment
READ(1,1001)keyname, agr, keycomment
READ(1,1001)keyname, drho, keycomment
READ(1,1001)keyname, dust2gas, keycomment
READ(1,1001)keyname, ebed, keycomment
READ(1,1000)keyname, tunneling, keycomment
READ(1,1000)keyname, barrier_tunneling, keycomment
READ(1,1001)keyname, barrier_tunneling_w, keycomment
READ(1,1000)keyname, btw_ch3oh_only
READ(1,1000)keyname, hop_act_competition, keycomment
READ(1,1000)keyname, SHINGLEDECKER_TUNN
READ(1,1000)keyname, eqtype, keycomment
READ(1,1001)keyname, sitedens, keycomment
READ(1,1000)keyname, bulk_chemistry
READ(1,1001)keyname, ebed_factor
READ(1,1001)keyname, bulk_diff_slowdown
READ(1,1001)keyname, effsurfmass
READ(1,*)
READ(1,1003)keyname, init_non_zero, keycomment
READ(1,*)
READ(1,1003)keyname, timesteps, keycomment
READ(1,1001)keyname, tstart, keycomment
READ(1,1001)keyname, tend, keycomment

ALLOCATE(timesteps_out(timesteps))

READ(1,*)
READ(1,1003)keyname, n_det_spec, keycomment
ALLOCATE(s_det_study(n_det_spec))

DO i = 1, n_det_spec
  READ(1,'(a10)')s_det_study(i)%name
ENDDO

CLOSE(1)

IF (delta_rho == 1) CALL read_d_profile
IF (delta_t == 1) CALL read_t_profile
IF (delta_tdust == 1) CALL read_tdust_profile
IF (delta_g0 == 1) CALL read_g0_profile
IF (delta_avst == 1) CALL read_avst_profile
IF (delta_avis == 1) CALL read_avis_profile
IF (delta_zetacr == 1) CALL read_zetacr_profile
IF (delta_zetax == 1) CALL read_zetax_profile
IF (delta_selfshield == 1) CALL read_selfshield_profile


! Added by C. N. Shingledecker
! 19.9.2018
Ngas  = 0
Nsurf = 0
Nbulk = 0
IF (init_non_zero .EQ. 1) THEN
  ! Count number of lines in the init_gas_ab.inp file
  io = 0
  OPEN(2,FILE='init_gas_ab.inp',STATUS='OLD',IOSTAT=io)
  DO
    READ(2,*,IOSTAT=io)
    IF (io .NE. 0) EXIT
    Ngas = Ngas + 1
  ENDDO
  CLOSE(2)

  ! Count number of lines in the init_surf_ab.inp file
  io = 0
  OPEN(3,FILE='init_surf_ab.inp',STATUS='OLD',IOSTAT=io)
  DO
    READ(3,*,IOSTAT=io)
    IF (io .NE. 0) EXIT
    Nsurf = Nsurf + 1
  ENDDO
  CLOSE(3)

  ! Count number of lines in the init_bulk_ab.inp file
  io = 0
  OPEN(4,FILE='init_bulk_ab.inp',STATUS='OLD',IOSTAT=io)
  DO
    READ(4,*,IOSTAT=io)
    IF (io .NE. 0) EXIT
    Nbulk = Nbulk + 1
  ENDDO
  CLOSE(4)

  init_non_zero = Ngas + Nsurf + Nbulk
ENDIF

ALLOCATE(s_init(init_non_zero))

! Added by C. N. Shingledecker
! 19.9.2018
IF (init_non_zero .GT. 1) THEN
  ! Read Gas Species File Again
  io = 0
  OPEN(2,FILE='init_gas_ab.inp',STATUS='OLD',IOSTAT=io)
  DO i = 1,Ngas
    READ(2,1004)s_init(i)%name, s_init(i)%abundance, keycomment
  ENDDO
  CLOSE(2)

  ! Read Surface Species File Again
  io = 0
  OPEN(3,FILE='init_surf_ab.inp',STATUS='OLD',IOSTAT=io)
  DO i = Ngas+1,Ngas+Nsurf
    READ(3,1004)s_init(i)%name, s_init(i)%abundance, keycomment
  ENDDO
  CLOSE(3)

  ! Read Bulk Species File Again
  io = 0
  OPEN(4,FILE='init_bulk_ab.inp',STATUS='OLD',IOSTAT=io)
  DO i = Ngas+Nsurf+1,init_non_zero
    READ(4,1004)s_init(i)%name, s_init(i)%abundance, keycomment
  ENDDO
  CLOSE(4)
ENDIF

!Check electric neutrality of the medium:

positive_ion_abundance = 0.0D0
negative_ion_and_electron_abundance = 0.0D0

DO i = 1, init_non_zero
    IF (s_init(i)%name(LEN_TRIM(s_init(i)%name):LEN_TRIM(s_init(i)%name)) == '+') positive_ion_abundance = positive_ion_abundance + s_init(i)%abundance
    IF (s_init(i)%name(LEN_TRIM(s_init(i)%name):LEN_TRIM(s_init(i)%name)) == '-') negative_ion_and_electron_abundance = negative_ion_and_electron_abundance + s_init(i)%abundance
ENDDO

IF (abs(positive_ion_abundance - negative_ion_and_electron_abundance) > epsilon(positive_ion_abundance - negative_ion_and_electron_abundance)) THEN
    PRINT*, positive_ion_abundance, negative_ion_and_electron_abundance, positive_ion_abundance - negative_ion_and_electron_abundance
    PRINT*, 'The medium is not electrically neutral! '
    PRINT*, 'Correcting election (e-) abundance...'

    electron_correction = 0.0D0

    DO i = 1, init_non_zero
        IF (s_init(i)%name(LEN_TRIM(s_init(i)%name):LEN_TRIM(s_init(i)%name)) == '+') electron_correction = electron_correction + s_init(i)%abundance
        IF (s_init(i)%name(LEN_TRIM(s_init(i)%name):LEN_TRIM(s_init(i)%name)) == '-') electron_correction = electron_correction - s_init(i)%abundance
        IF (s_init(i)%name == 'e-') i_e = i
    ENDDO

    s_init(i_e)%abundance = s_init(i_e)%abundance + electron_correction

ENDIF


RETURN
!Error statements:
100 PRINT*, 'File model.inp not found!'
STOP
!Format statements:
1000 FORMAT(a23, 1x, i1, 10x, a40)
1001 FORMAT(a23, 1x, e10.8, 1x, a40)
1002 FORMAT(a23, 1x, a80)
1003 FORMAT(a23, 1x, i3, 8x ,a40)
1004 FORMAT(a10, 14x, e14.12, 1x, a40)
1007 FORMAT(a23, 1x, l1, 10x, a40)

END SUBROUTINE read_model_setup

SUBROUTINE read_d_profile
IMPLICIT NONE

n_d_steps = 0
OPEN(1,FILE='d_profile',STATUS='OLD',ERR=102)

DO WHILE (.NOT. EOF(1))
  n_d_steps = n_d_steps + 1
  READ(1,'(1pe12.4,1x,1pe12.4)')time_d_array(n_d_steps), d_array(n_d_steps)
END DO

WRITE(*,'(a35,1pe12.4,a11,1pe12.4,a6)') 'read_d_profile: initial density = ',d_array(1),' at time = ',time_d_array(1),' years'
WRITE(*,'(a35,1pe12.4,a11,1pe12.4,a6)') 'read_d_profile: final density   = ',d_array(n_d_steps),' at time = ',time_d_array(n_d_steps),' years'
PRINT*, 'read_d_profile: amount of steps = ',n_d_steps

CLOSE(1)

RETURN
102 PRINT*, 'File d_profile not found!'

END SUBROUTINE read_d_profile

SUBROUTINE read_t_profile
IMPLICIT NONE

n_t_steps = 0

OPEN(1,FILE='t_profile',STATUS='OLD',ERR=103)

DO WHILE (.NOT. EOF(1))
  n_t_steps = n_t_steps + 1
  READ(1,'(1pe12.4,1x,1pe12.4)')time_t_array(n_t_steps), t_array(n_t_steps)
END DO

WRITE(*,'(a39,1pe12.4,a11,1pe12.4,a6)') 'read_t_profile: initial temperature = ',t_array(1),' at time = ',time_t_array(1),' years'
WRITE(*,'(a39,1pe12.4,a11,1pe12.4,a6)') 'read_t_profile: final temperature   = ',t_array(n_t_steps),' at time = ',time_t_array(n_t_steps),' years'
PRINT*, 'read_t_profile: amount of steps = ',n_t_steps

CLOSE(1)

RETURN
103 PRINT*, 'File t_profile not found!'

END SUBROUTINE read_t_profile

SUBROUTINE read_tdust_profile
IMPLICIT NONE

n_tdust_steps = 0

OPEN(1,FILE='tdust_profile',STATUS='OLD',ERR=1033)

DO WHILE (.NOT. EOF(1))
  n_tdust_steps = n_tdust_steps + 1
  READ(1,'(1pe12.4,1x,1pe12.4)')time_tdust_array(n_tdust_steps), tdust_array(n_tdust_steps)
END DO

WRITE(*,'(a39,1pe12.4,a11,1pe12.4,a6)') 'read_tdust_profile: initial temperature = ',tdust_array(1),' at time = ',time_tdust_array(1),' years'
WRITE(*,'(a39,1pe12.4,a11,1pe12.4,a6)') 'read_tdust_profile: final temperature   = ',tdust_array(n_tdust_steps),' at time = ',time_tdust_array(n_tdust_steps),' years'
PRINT*, 'read_tdust_profile: amount of steps = ',n_tdust_steps

CLOSE(1)

RETURN
1033 PRINT*, 'File tdust_profile not found!'

END SUBROUTINE read_tdust_profile

SUBROUTINE read_g0_profile
IMPLICIT NONE

n_g0_steps = 0

OPEN(1,FILE='g0_profile',STATUS='OLD',ERR=104)

DO WHILE (.NOT. EOF(1))
  n_g0_steps = n_g0_steps + 1
  READ(1,'(1pe12.4,1x,1pe12.4)')time_g0_array(n_g0_steps), g0_array(n_g0_steps)
END DO

WRITE(*,'(a39,1pe12.4,a11,1pe12.4,a6)') 'read_g0_profile: initial G0 = ',g0_array(1),' at time = ',time_g0_array(1),' years'
WRITE(*,'(a39,1pe12.4,a11,1pe12.4,a6)') 'read_g0_profile: final G0 = ',g0_array(n_g0_steps),' at time = ',time_g0_array(n_g0_steps),' years'
PRINT*, 'read_g0_profile: amount of steps = ',n_g0_steps

CLOSE(1)

RETURN
104 PRINT*, 'File g0_profile not found!'

END SUBROUTINE read_g0_profile

SUBROUTINE read_avst_profile
IMPLICIT NONE

n_avst_steps = 0

OPEN(1,FILE='avst_profile',STATUS='OLD',ERR=105)

DO WHILE (.NOT. EOF(1))
  n_avst_steps = n_avst_steps + 1
  READ(1,'(1pe12.4,1x,1pe12.4)')time_avst_array(n_avst_steps), avst_array(n_avst_steps)
END DO

WRITE(*,'(a39,1pe12.4,a11,1pe12.4,a6)') 'read_avst_profile: initial AvSt = ',avst_array(1),' at time = ',time_avst_array(1),' years'
WRITE(*,'(a39,1pe12.4,a11,1pe12.4,a6)') 'read_avst_profile: final AvSt = ',avst_array(n_avst_steps),' at time = ',time_avst_array(n_avst_steps),' years'
PRINT*, 'read_avst_profile: amount of steps = ',n_avst_steps

CLOSE(1)

RETURN
105 PRINT*, 'File avst_profile not found!'

END SUBROUTINE read_avst_profile

SUBROUTINE read_avis_profile
IMPLICIT NONE

n_avis_steps = 0

OPEN(1,FILE='avis_profile',STATUS='OLD',ERR=106)

DO WHILE (.NOT. EOF(1))
  n_avis_steps = n_avis_steps + 1
  READ(1,'(1pe12.4,1x,1pe12.4)')time_avis_array(n_avis_steps), avis_array(n_avis_steps)
END DO

WRITE(*,'(a39,1pe12.4,a11,1pe12.4,a6)') 'read_avis_profile: initial AvIS = ',avis_array(1),' at time = ',time_avis_array(1),' years'
WRITE(*,'(a39,1pe12.4,a11,1pe12.4,a6)') 'read_avis_profile: final AvIS = ',avis_array(n_avis_steps),' at time = ',time_avis_array(n_avis_steps),' years'
PRINT*, 'read_avis_profile: amount of steps = ',n_avis_steps

CLOSE(1)

RETURN
106 PRINT*, 'File avis_profile not found!'

END SUBROUTINE read_avis_profile

SUBROUTINE read_zetacr_profile
IMPLICIT NONE

n_zetacr_steps = 0

OPEN(1,FILE='zetacr_profile',STATUS='OLD',ERR=107)

DO WHILE (.NOT. EOF(1))
  n_zetacr_steps = n_zetacr_steps + 1
  READ(1,'(1pe12.4,1x,1pe12.4)')time_zetacr_array(n_zetacr_steps), zetacr_array(n_zetacr_steps)
END DO

WRITE(*,'(a39,1pe12.4,a11,1pe12.4,a6)') 'read_zetacr_profile: initial ZetaCR = ',zetacr_array(1),' at time = ',time_zetacr_array(1),' years'
WRITE(*,'(a39,1pe12.4,a11,1pe12.4,a6)') 'read_zetacr_profile: final ZetaCR = ',zetacr_array(n_zetacr_steps),' at time = ',time_zetacr_array(n_zetacr_steps),' years'
PRINT*, 'read_zetacr_profile: amount of steps = ',n_zetacr_steps

CLOSE(1)

RETURN
107 PRINT*, 'File zetacr_profile not found!'

END SUBROUTINE read_zetacr_profile

SUBROUTINE read_zetax_profile
IMPLICIT NONE

n_zetax_steps = 0

OPEN(1,FILE='zetax_profile',STATUS='OLD',ERR=108)

DO WHILE (.NOT. EOF(1))
  n_zetax_steps = n_zetax_steps + 1
  READ(1,'(1pe12.4,1x,1pe12.4)')time_zetax_array(n_zetax_steps), zetax_array(n_zetax_steps)
END DO

WRITE(*,'(a39,1pe12.4,a11,1pe12.4,a6)') 'read_zetax_profile: initial zetax = ',zetax_array(1),' at time = ',time_zetax_array(1),' years'
WRITE(*,'(a39,1pe12.4,a11,1pe12.4,a6)') 'read_zetax_profile: final zetax = ',zetax_array(n_zetax_steps),' at time = ',time_zetax_array(n_zetax_steps),' years'
PRINT*, 'read_zetax_profile: amount of steps = ',n_zetax_steps

CLOSE(1)

RETURN
108 PRINT*, 'File zetax_profile not found!'

END SUBROUTINE read_zetax_profile

SUBROUTINE read_selfshield_profile
IMPLICIT NONE

n_selfshield_steps = 0

OPEN(1,FILE='selfshield_profile',STATUS='OLD',ERR=109)

DO WHILE (.NOT. EOF(1))
  n_selfshield_steps = n_selfshield_steps + 1
  READ(1,'(1pe12.4,4(1x,1pe12.4))')time_selfshield_array(n_selfshield_steps), fh2is_array(n_selfshield_steps), fcois_array(n_selfshield_steps), fh2st_array(n_selfshield_steps), fcost_array(n_selfshield_steps)
END DO

PRINT*, 'read_selfshield_profile: amount of steps = ',n_selfshield_steps

CLOSE(1)

RETURN
109 PRINT*, 'File selfshield_profile not found!'

END SUBROUTINE read_selfshield_profile

END
