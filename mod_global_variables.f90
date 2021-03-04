  
MODULE global_variables
IMPLICIT NONE

INTEGER, PARAMETER :: wp = KIND(1.0d0)

TYPE species
  CHARACTER*10                        :: name
  REAL(wp)                            :: weight, edes, abundance, racc, rdes, enthalpia, bdiffrate
  REAL(wp)                            :: frac_abundance
  REAL(wp), DIMENSION(:), ALLOCATABLE :: abundance_out
  INTEGER                             :: idx, gas_idx, enthalpia_known, natoms
!  INTEGER                             :: nH,nC,nN,nO,nS,nP,nD,nF,nHe,nFe,nSi,nNa,nMg,nCl
!  INTEGER                             :: nPlusCharge,nMinusCharge
END TYPE species

TYPE reaction
  INTEGER      :: idx
  INTEGER      :: rtype, exothermicity_known
  CHARACTER*10 :: r1, r2, p1, p2, p3, p4, p5
  INTEGER      :: ir1, ir2, ip1, ip2, ip3, ip4, ip5
  REAL(wp)     :: alpha, beta, gamma
  REAL(wp)     :: rate, exothermicity
END TYPE reaction

!Physical constants:
REAL(wp), PARAMETER                          :: aMp=1.6724E-24, aMe=9.1095E-28, amu=1.43d0, hp=6.6262E-27, ak_B=1.3807E-16, year=3.155E+07, PI=3.141593d0, albedo_UV=0.5d0, yield=1.0E-5
INTEGER, SAVE                                :: sing_mult, &
                                                is_disk_model, &
                                                radiolysis, &
                                                suprathermal, &
                                                tunneling, &
                                                barrier_tunneling, &
                                                eqtype, bulk_chemistry, &
                                                init_non_zero, &
                                                first_surfreact, &
                                                first_bulkreact, &
                                                first_surf_spec, &
                                                n_surf_spec, &
                                                n_surf_react, &
                                                is_pr_precalculated, &
                                                hop_act_competition, &
                                                btw_ch3oh_only, &
                                                shingledecker_tunn, &
                                                model_experiment,&
                                                fast_bulk,&
                                                fast_atoms,&
                                                photoexc,&
                                                photoion,&
                                                FIXED_DVAL,&
                                                FIXED_NU,&
                                                DISABLE_DESORB
INTEGER, SAVE                                :: delta_rho, delta_t, delta_tdust, delta_g0, delta_avst, delta_avis, delta_zetacr, delta_zetax, delta_selfshield
INTEGER, SAVE                                :: n_d_steps, n_t_steps, n_tdust_steps, n_G0_steps, n_avst_steps, n_avis_steps, n_zetacr_steps, n_zetax_steps, n_selfshield_steps
REAL(wp), SAVE                               :: rho, rho_ice, ice_thick, t, tdust, G0_stellar, AvSt, AvIS, ZetaCR, ZetaX, agr, drho, dust2gas, ebed, ebed_factor, sitedens, nsites, nml, nml_old, nml_max, n_s_ml, bulk_diff_slowdown, barrier_tunneling_w, ph_yield, phi_exp,se_exp,DVAL,EXTFAC
REAL(wp), SAVE                               :: effsurfmass, gdens, ddens, tstart, tend, told, rtol, atol
REAL(wp), DIMENSION(2000)                    :: d_array, t_array, tdust_array, G0_array, avst_array, avis_array, zetacr_array, zetax_array, fh2is_array, fcois_array, fh2st_array, fcost_array
REAL(wp), DIMENSION(2000)                    :: time_d_array, time_t_array, time_tdust_array, time_G0_array, time_avst_array, time_avis_array, time_zetacr_array, time_zetax_array, time_selfshield_array
CHARACTER*80, SAVE                           :: chem_file, outfile_nameprefix, photorate_nameprefix
INTEGER, SAVE                                :: nspecies, nreactions, timesteps, n_det_spec
INTEGER, SAVE                                :: nr,nz, curpoint, fdp, ldp, rde, myunit
INTEGER, SAVE                                :: des_reactive_type
INTEGER                                      :: des_t, des_crp, des_photon
REAL(wp), SAVE                               :: des_reactive
REAL(wp), SAVE                               :: rs, zs, fh2_is, fh2_st, fco_is, fco_st
REAL(wp), SAVE                               :: dtran, smol, bmol
REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: reaction_importance !species, reaction, time_moment
REAL(wp), DIMENSION(:,:), ALLOCATABLE        :: abundances_bulk
REAL(wp), DIMENSION(:), ALLOCATABLE          :: timesteps_out, timesteps_nml, mre_terms
REAL(wp), DIMENSION(:,:), ALLOCATABLE        :: rd_v2_terms
TYPE (species), DIMENSION(:), ALLOCATABLE    :: s, s_init, s_det_study
TYPE (reaction), DIMENSION(:), ALLOCATABLE   :: r
CHARACTER(LEN=4), DIMENSION(9)              :: radical_names = (/"bO  ","bOH ","bHO2","bHO3","bH  ","bHS ","bNS ","bHSO","bCS "/)
REAL(wp)                                     :: initial_water, trial_nu

END