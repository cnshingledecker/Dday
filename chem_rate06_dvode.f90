PROGRAM chem_rate06_dvode
USE global_variables
USE read_model
USE calculate_rates
USE read_rate06
USE run_dvode
USE save_results

IMPLICIT NONE

INTEGER :: i

rde = 0
OPEN(newunit=myunit,file='rd_eff.txt')

CALL read_model_setup
CALL read_rate06database

IF (MODEL_EXPERIMENT.EQ.0) THEN
  nsites = sitedens*2.0d0*pi*agr**2.0d0
ELSE
  nsites = sitedens*ICE_AREA ! number of sites on 1um^2
ENDIF

ALLOCATE (reaction_importance(n_det_spec,nreactions,timesteps))
call save_Results_bulk

print*, 'Sing_Mult           = ', sing_mult
print*, 'is_disk_model       = ', is_disk_model
print*, 'radiolysis          = ', radiolysis
print*, 'chem_file           = ', chem_file(1:LEN_TRIM(chem_file))
print*, 'n_s_ml              = ', n_s_ml
print*, 'Density             = ', rho
print*, 'Temperature         = ', t
print*, 'G0_stellar          = ', G0_stellar
print*, 'AvSt                = ', AvSt
print*, 'AvIS                = ', AvIS
print*, 'ZetaCR              = ', ZetaCR
print*, 'ZetaX               = ', ZetaX
print*, 'des_t               = ', des_t
print*, 'des_crp             = ', des_crp
print*, 'des_photon          = ', des_photon
print*, 'ph_yield            = ', ph_yield
print*, 'des_reactive        = ', des_reactive
print*, 'des_reactive_type   = ', des_reactive_type
print*, 'agr                 = ', agr
print*, 'drho                = ', drho
print*, 'dust2gas            = ', dust2gas
print*, 'ebed                = ', ebed
print*, 'tunneling           = ', tunneling
print*, 'barrier_tunneling   = ', barrier_tunneling
print*, 'barrier_tunneling_w = ', barrier_tunneling_w
print*, 'btw_ch3oh_only      = ', btw_ch3oh_only
print*, 'hop_act_competition = ', hop_act_competition
print*, 'eqtype              = ', eqtype
print*, 'sitedens            = ', sitedens
print*, 'ini_non_zero        = ', init_non_zero
print*, 'nspecies            = ', nspecies
print*, 'nreactions          = ', nreactions

!pause

DO i = 1, init_non_zero
  print*, s_init(i)%name, s_init(i)%abundance
ENDDO

CALL calc_rates(0.d0)

SELECT CASE (eqtype)
  CASE (1)
    CALL run_dvode_solver_sparse
  CASE (2)
    CALL run_dvode_solver_sparse_mre
  CASE (3)
    PRINT*, 'Moment equations are not implemented yet!'
  CASE (4)
      CALL run_dvode_solver
  CASE DEFAULT
    PRINT*, 'Unknown type of equations: ', eqtype
END SELECT


CALL save_results_shingledecker

DEALLOCATE (s)
DEALLOCATE (r)
DEALLOCATE (s_init)
DEALLOCATE (timesteps_out)
DEALLOCATE (s_det_study)
DEALLOCATE (mre_terms)
DEALLOCATE (rd_v2_terms)
DEALLOCATE (reaction_importance)

END
