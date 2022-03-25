MODULE calculate_rates
USE global_variables

IMPLICIT NONE

CONTAINS

SUBROUTINE calc_rates(tcur)
IMPLICIT NONE
REAL(wp) Stick, Stick_H, Cion, Vth, anu0, anu1, ak_photo, ak_therm, ak_crp, tcrp, Rdiff0, Rdiff1, akbar, amur, dust2gas_corr, tcur
REAL(wp) alpha,beta,gamma,T0,tunn
INTEGER i
CHARACTER(LEN=10)         :: reactant_names(2)
CHARACTER(LEN=10)         :: r1,r2
REAL(wp)                  :: f1,f2,rate,r_alpha


IF (delta_t==1) T = get_parameter(tcur/year, n_t_steps, time_t_array, t_array)
IF (delta_tdust==1) Tdust = get_parameter(tcur/year, n_tdust_steps, time_tdust_array, tdust_array)
IF (delta_g0==1) G0_stellar = get_parameter(tcur/year, n_g0_steps, time_g0_array, g0_array)
IF (delta_avst==1) AvSt = get_parameter(tcur/year, n_avst_steps, time_avst_array, avst_array)
IF (delta_avis==1) AvIS = get_parameter(tcur/year, n_avis_steps, time_avis_array, avis_array)
IF (delta_zetacr==1) ZetaCR = get_parameter(tcur/year, n_zetacr_steps, time_zetacr_array, zetacr_array)
IF (delta_zetax==1) ZetaX = get_parameter(tcur/year, n_zetax_steps, time_zetax_array, zetax_array)
IF (delta_selfshield==1) THEN
    fh2_is = get_parameter(tcur/year, n_selfshield_steps, time_selfshield_array, fh2is_array)
    fco_is = get_parameter(tcur/year, n_selfshield_steps, time_selfshield_array, fcois_array)
    fh2_st = get_parameter(tcur/year, n_selfshield_steps, time_selfshield_array, fh2st_array)
    fco_st = get_parameter(tcur/year, n_selfshield_steps, time_selfshield_array, fcost_array)
ENDIF

dust2gas_corr = max(1.d-10,dust2gas)
IF (delta_rho==1) gdens = get_parameter(tcur/year, n_d_steps, time_d_array, d_array)
ddens = gdens * dust2gas_corr / (4.0d0/3.0d0*PI*agr*agr*agr*drho) * aMp * amu
Stick = 1.0D0
Stick_H = 1.0D0 / (1D0+4D-2*DSQRT(Tdust+Tdust)+2D-3*Tdust+8D-6*Tdust*Tdust)
Cion = 1.0D0 + 1.671D-3 / agr / T

CALL shield(T,AvSt,fh2_st,fco_st)
CALL shield(T,AvIS,fh2_is,fco_is)

ak_photo = 0.d0
ak_therm = 0.d0
ak_crp   = 0.d0

DO i = 1, nreactions
    SELECT CASE (r(i)%rtype)
        CASE (1) !Two-body molecule(ion) - molecule(ion) gas-phase reaction
            r(i)%rate = r(i)%alpha*(T/300.0d0)**r(i)%beta*dexp(-r(i)%gamma/T) * ddens
        CASE (2) !Cosmic ray ionization reaction (CRP)
            IF ( MODEL_EXPERIMENT .EQ. 0 ) THEN
              r(i)%rate = r(i)%alpha/1.3E-17*(zetaCR+zetaX)
            ELSE
              r(i)%rate = 0.0d0
            ENDIF
        CASE (3) !Photoionization reaction (PHOTON)
            IF ( MODEL_EXPERIMENT .EQ. 0 ) THEN
              r(i)%rate = r(i)%alpha*(DEXP(-r(i)%gamma*AvSt)*G0_stellar+DEXP(-r(i)%gamma*AvIS))
            ELSE
              r(i)%rate = 0.0d0
            ENDIF
        CASE (4) !Photoionization, H2 self-shielded
            IF ( MODEL_EXPERIMENT .EQ. 0 ) THEN
              r(i)%rate = r(i)%alpha*(G0_Stellar*fh2_St+fh2_IS)
            ELSE
              r(i)%rate = 0.0d0
            ENDIF
        CASE (5) !Photoionization, CO self-shielded
            IF ( MODEL_EXPERIMENT .EQ. 0 ) THEN
              r(i)%rate = r(i)%alpha*(G0_Stellar*fco_St+fco_IS)
            ELSE
              r(i)%rate = 0.0d0
            ENDIF
        CASE (6) !Cosmic ray-induced photoreaction (CRPHOT)
            IF ( MODEL_EXPERIMENT .EQ. 0 ) THEN
              r(i)%rate = r(i)%alpha/1.3E-17*(zetaCR+zetaX)*r(i)%gamma/(1.0d0-albedo_UV)
            ELSE
              r(i)%rate = 0.0d0
            ENDIF
        CASE (7) !Ions plus negatively charged grains
            Vth = DSQRT(8.0d0*ak_B*T/PI/s(r(i)%ir1)%weight/aMp)
            r(i)%rate = r(i)%alpha*PI*agr*agr*Vth*Stick*Cion*ddens
        CASE (8) !Ions plus neutral grains
            Vth = DSQRT(8.0d0*ak_B*T/PI/s(r(i)%ir1)%weight/aMp)
            r(i)%rate = r(i)%alpha*PI*agr*agr*Vth*Stick*ddens
        CASE (9) !Collisions of the electrons with the neutral grains
            Vth = DSQRT(8.0D0*ak_B*T/PI/aMe)
            r(i)%rate = PI*agr*agr*Vth*1.329D0*DEXP(-T/20.0D0)*ddens
        CASE (10) !Collisions of the electrons with the positively charged grains
            Vth = DSQRT(8.0D0*ak_B*T/PI/aMe)
            r(i)%rate = PI*agr*agr*Vth*1.329D0*DEXP(-T/20.0D0)*Cion*ddens
        CASE (11) !Accretion (FREEZE)
            Vth = DSQRT(8.0d0*ak_B*T/PI/s(r(i)%ir1)%weight/aMp)
            r(i)%rate = r(i)%alpha*PI*agr*agr*Vth*Stick*ddens
            IF (r(i)%r1=='H') r(i)%rate = r(i)%alpha*PI*agr*agr*Vth*Stick_H*ddens
            s(r(i)%ir1)%racc = r(i)%rate
        CASE (12) !Desorption (DESORB)
            anu0 = DSQRT(2.d0*sitedens*ak_B/aMp*r(i)%gamma/PI/PI/s(r(i)%ir1)%weight)
            ! Photodesorption:
            IF ( des_photon==1 .AND. MODEL_EXPERIMENT .EQ. 0 ) THEN
              ak_photo = (ph_yield+0.13d0*DEXP(-336.0d0/Tdust))*3.38e-8*(G0_stellar*dexp(-2.0d0*AvSt)+dexp(-2.0d0*AvIS))
            ELSEIF ( des_photon==1 .AND. MODEL_EXPERIMENT .EQ. 1 ) THEN
              ak_photo = ph_yield*PHI_EXP
            ENDIF
            ! Thermal desorption:
            IF (des_t==1) ak_therm = anu0*DEXP(-r(i)%gamma/Tdust)
            ! CRP desorption:
            ! Le'ger et al. 1985 (CRP heating of 1000A grain, Eq. 1a):
            IF (des_crp==1) ak_crp = 3.16d-19*anu0*DEXP(-r(i)%gamma/70.d0)*(zetaCR/1.3d-17) !2.431E-2*zetaCR*anu0*DEXP(-r(i)%gamma/tcrp)
            IF ( DISABLE_DESORB .EQ. 0 ) THEN
              r(i)%rate = r(i)%alpha*(ak_photo + ak_therm + ak_crp)
            ELSE
              r(i)%rate = 0.0d0
            ENDIF
            s(r(i)%ir1)%rdes = r(i)%rate
        CASE (13) !Surface two-body reaction
            ! Compute their vibrational frequencies:
            anu0 = DSQRT(sitedens*ak_B/aMp*s(r(i)%ir1)%edes/PI/PI/s(r(i)%ir1)%weight)
            anu1 = DSQRT(sitedens*ak_B/aMp*s(r(i)%ir2)%edes/PI/PI/s(r(i)%ir2)%weight)
            ! Compute diffusion rates of the reactants (thermal hopping and quantum tunneling for the model of HHL 1992):
            Rdiff0 = anu0/(4.0d0*pi*agr**2.0d0*sitedens)*DEXP(-ebed*s(r(i)%ir1)%edes/Tdust)
            Rdiff1 = anu1/(4.0d0*pi*agr**2.0d0*sitedens)*DEXP(-ebed*s(r(i)%ir2)%edes/Tdust)


            !Diffusion barriers:

            IF (tunneling==1) THEN
                IF (s(r(i)%ir1)%name=='gH' .OR. s(r(i)%ir1)%name=='gH2') THEN
                    amur = s(r(i)%ir1)%weight*aMp
                    akbar=DEXP(-4D0*PI/hp*1D-8*DSQRT(2D0*amur*ebed*s(r(i)%ir1)%edes*ak_B))
                    Rdiff0 = anu0/(2.0d0*pi*agr**2.0d0*sitedens)*akbar
                ENDIF
                IF (s(r(i)%ir2)%name=='gH' .OR. s(r(i)%ir2)%name=='gH2') THEN
                    amur = s(r(i)%ir2)%weight*aMp
                    akbar=DEXP(-4D0*PI/hp*1D-8*DSQRT(2D0*amur*ebed*s(r(i)%ir2)%edes*ak_B))
                    Rdiff1 = anu1/(2.0d0*pi*agr**2.0d0*sitedens)*akbar
                ENDIF
                ! ED = 1040 K = 520 K / 0.5, with 520 K being the best fit barrier height
                ! a = 0.7 Angstrom
                ! See Minissale et al. 2014
                IF (s(r(i)%ir1)%name=='gO') THEN
                  amur   = s(r(i)%ir2)%weight*aMp
                  akbar  = DEXP(-4D0*PI/hp*0.7D-8*DSQRT(2D0*amur*ebed*1040.0*ak_B))
                  Rdiff0 = anu0/(2.0d0*pi*agr**2.0d0*sitedens)*akbar
                ELSE IF (s(r(i)%ir2)%name=='gO') THEN
                  amur   = s(r(i)%ir2)%weight*aMp
                  akbar  = DEXP(-4D0*PI/hp*0.7D-8*DSQRT(2D0*amur*ebed*1040.0*ak_B))
                  Rdiff1 = anu1/(2.0d0*pi*agr**2.0d0*sitedens)*akbar
                ENDIF
            ENDIF


            !Reaction activation barriers:
            akbar= DEXP(-r(i)%gamma/Tdust)

            IF (barrier_tunneling==1) THEN
                IF (s(r(i)%ir1)%name=='gH' .OR. s(r(i)%ir1)%name=='gH2' .OR. s(r(i)%ir2)%name=='gH' .OR. s(r(i)%ir2)%name=='gH2') THEN
                    amur = s(r(i)%ir1)%weight*s(r(i)%ir2)%weight/(s(r(i)%ir1)%weight+s(r(i)%ir2)%weight)*aMp
                    IF (btw_ch3oh_only==1) THEN
                        akbar=DEXP(-4d0*PI/hp*1.0D-8*DSQRT(2d0*amur*r(i)%gamma*ak_B))
                        IF (s(r(i)%ir1)%name=='gH' .AND. s(r(i)%ir2)%name=='gCO') akbar=DEXP(-4d0*PI/hp*barrier_tunneling_w*DSQRT(2d0*amur*r(i)%gamma*ak_B))
                        IF (s(r(i)%ir1)%name=='gH' .AND. s(r(i)%ir2)%name=='gH2CO') akbar=DEXP(-4d0*PI/hp*barrier_tunneling_w*DSQRT(2d0*amur*r(i)%gamma*ak_B))
                    ELSE
                        akbar=DEXP(-4d0*PI/hp*barrier_tunneling_w*DSQRT(2d0*amur*r(i)%gamma*ak_B))
                    ENDIF
                ENDIF
            ENDIF

            IF (hop_act_competition==1 .AND. shingledecker_tunn .EQ. 0) THEN
                IF (r(i)%gamma/=0.d0) THEN
                  akbar = (akbar*anu0+akbar*anu1)/(akbar*anu0+akbar*anu1+Rdiff0*(2.0d0*pi*agr**2.0d0*sitedens)+Rdiff1*(2.0d0*pi*agr**2.0d0*sitedens))
!                  IF (s(r(i)%ir1)%name=='gH' .OR. s(r(i)%ir1)%name=='gH2' .OR. s(r(i)%ir2)%name=='gH' .OR. s(r(i)%ir2)%name=='gH2') THEN
!                    PRINT *, r(i)%r1," = H or H2, and/or, ",r(i)%r2, " with akbar= ",akbar
!                  ELSE
!                    PRINT *, r(i)%r1," NOT EQUAL TO H or H2, and/or, ",r(i)%r2, " with akbar= ",akbar
!                  ENDIF
                ENDIF
            ELSEIF ( hop_act_competition.EQ. 1 .AND. shingledecker_tunn .EQ. 1 ) THEN
              ! IF one of the reactants is H or H2, the reaction probability is (tunnel rate + thermal rate)/(tunnel rate + thermal rate + diffrate1 + diffrate2)
              IF (s(r(i)%ir1)%name=='gH' .OR. s(r(i)%ir1)%name=='gH2' .OR. s(r(i)%ir2)%name=='gH' .OR. s(r(i)%ir2)%name=='gH2') THEN
                IF (r(i)%gamma/=0.d0) THEN
                  akbar = (&
                                               anu0*DEXP(-r(i)%gamma/Tdust)+&
                                               anu1*DEXP(-r(i)%gamma/Tdust)+&
                                               akbar*anu0+&
                                               akbar*anu1&
                                              )/&
                                              (&
                                               anu0*DEXP(-r(i)%gamma/Tdust)+&
                                               anu1*DEXP(-r(i)%gamma/Tdust)+&
                                               akbar*anu0+&
                                               akbar*anu1+&
                                               Rdiff0*(2.0d0*pi*agr**2.0d0*sitedens)+&
                                               Rdiff1*(2.0d0*pi*agr**2.0d0*sitedens)&
                                              )
!                  PRINT *, r(i)%r1," = H or H2, and/or, ",r(i)%r2, " with akbar= ",akbar
                ENDIF
              ! If there is no tunnel rate, the only thermal rate is used in the numerator
              ELSE
                IF (r(i)%gamma/=0.d0) THEN
                  akbar = (akbar*anu0+akbar*anu1)/(akbar*anu0+akbar*anu1+Rdiff0*(2.0d0*pi*agr**2.0d0*sitedens)+Rdiff1*(2.0d0*pi*agr**2.0d0*sitedens))
!                  PRINT *, r(i)%r1," NOT EQUAL TO H or H2, and/or, ",r(i)%r2, " with akbar= ",akbar
                ENDIF
              ENDIF
            ENDIF

            IF (shingledecker_tunn==1) THEN
              ! Special case for H2 + OH -> H2O + H
              IF ((s(r(i)%ir1)%name=='gH2' .AND. s(r(i)%ir2)%name=='gOH') .OR. &
                (s(r(i)%ir2)%name=='gH2' .AND. s(r(i)%ir1)%name=='gOH')) THEN
                alpha = 7.64e10
                beta  = 1.00
                gamma = 1339.9
                T0    = 153.2
                tunn  = alpha*((Tdust/300.0)**beta)*DEXP(-1*gamma*(Tdust + T0)/((Tdust**2) + (T0**2)))
                akbar= DEXP(-r(i)%gamma/Tdust)
                IF ( hop_act_competition .EQ. 1 ) THEN
                  ! Revised competition formula
                  akbar = (tunn+akbar*anu0+akbar*anu1)/(tunn+akbar*anu0+akbar*anu1+Rdiff0*&
                    (4.0d0*pi*agr**2.0d0*sitedens)+Rdiff1*(2.0d0*pi*agr**2.0d0*sitedens))
                ELSE
                  ! Just scale exponential
                  akbar = (tunn + (anu0+anu1)*DEXP(-r(i)%gamma/Tdust))/(anu0+anu1)
                  IF ( akbar .GT. 1.0 ) akbar = 1.0
                ENDIF
              ! Special case for H2 + O3 -> H2O + O2
!              ELSE IF ((s(r(i)%ir1)%name=='gH2' .AND. s(r(i)%ir2)%name=='gO3') .OR. &
!                (s(r(i)%ir2)%name=='gH2' .AND. s(r(i)%ir1)%name=='gO3')) THEN
!                alpha = 7.64e10
!                beta  = 1.00
!                gamma = 1339.9
!                T0    = 153.2
!                tunn  = alpha*((Tdust/300.0)**beta)*DEXP(-1*gamma*(Tdust + T0)/((Tdust**2) + (T0**2)))
!                amur  = s(r(i)%ir1)%weight*s(r(i)%ir2)%weight/(s(r(i)%ir1)%weight+s(r(i)%ir2)%weight)*aMp
!                akbar =DEXP(-4d0*PI/hp*barrier_tunneling_w*DSQRT(2d0*amur*r(i)%gamma*ak_B))
!                akbar = (tunn+akbar*anu0+akbar*anu1)/(tunn+akbar*anu0+akbar*anu1+Rdiff0*&
!                  (2.0d0*pi*agr**2.0d0*sitedens)+Rdiff1*(2.0d0*pi*agr**2.0d0*sitedens))
              ! Special case for H + H2O2 -> H2O + OH
              ELSE IF ((s(r(i)%ir1)%name=='gH2O2' .AND. s(r(i)%ir2)%name=='gH') .OR. &
                (s(r(i)%ir2)%name=='gH2O2' .AND. s(r(i)%ir1)%name=='gH')) THEN
                alpha = 1.51e10
                beta  = 0.86
                gamma = 1750.0
                T0    = 180.0
                tunn  = alpha*((Tdust/300.0)**beta)*DEXP(-1*gamma*(Tdust + T0)/((Tdust**2) + (T0**2)))
                akbar= DEXP(-r(i)%gamma/Tdust)
                IF ( hop_act_competition .EQ. 1 ) THEN
                  ! Revised competition formula
                  akbar = (tunn+akbar*anu0+akbar*anu1)/(tunn+akbar*anu0+akbar*anu1+Rdiff0*&
                    (4.0d0*pi*agr**2.0d0*sitedens)+Rdiff1*(2.0d0*pi*agr**2.0d0*sitedens))
                ELSE
                  ! Just scale exponential
                  akbar = (tunn + (anu0+anu1)*DEXP(-r(i)%gamma/Tdust))/(anu0+anu1)
                  IF ( akbar .GT. 1.0 ) akbar = 1.0
                ENDIF
              ! Special case for H + H2S ->  HS + H2
              ELSE IF ((s(r(i)%ir1)%name=='gH' .AND. s(r(i)%ir2)%name=='gH2S') .OR. &
                (s(r(i)%ir2)%name=='gH' .AND. s(r(i)%ir1)%name=='gH2S')) THEN
                alpha = 2.20111E11
                beta  = 0.475667
                gamma = 1399.51
                T0    = 180
                tunn  = alpha*((Tdust/300.0)**beta)*DEXP(-1*gamma*(Tdust + T0)/((Tdust**2) + (T0**2)))
                akbar= DEXP(-r(i)%gamma/Tdust)
                IF ( hop_act_competition .EQ. 1 ) THEN
                  ! Revised competition formula
                  akbar = (tunn+akbar*anu0+akbar*anu1)/(tunn+akbar*anu0+akbar*anu1+Rdiff0*&
                    (4.0d0*pi*agr**2.0d0*sitedens)+Rdiff1*(2.0d0*pi*agr**2.0d0*sitedens))
                ELSE
                  ! Just scale exponential
                  akbar = (tunn + (anu0+anu1)*DEXP(-r(i)%gamma/Tdust))/(anu0+anu1)
                  IF ( akbar .GT. 1.0 ) akbar = 1.0
                ENDIF
              ! Special case for H + CS -> HCS
              ELSE IF ((s(r(i)%ir1)%name=='gH' .AND. s(r(i)%ir1)%name=='gCS') .OR. &
                (s(r(i)%ir2)%name=='gH' .AND. s(r(i)%ir2)%name=='gCS')) THEN
                alpha = 2.05296E12
                beta  = -8.10946
                gamma = 984.986
                T0    = 50.1934
                tunn  = alpha*((Tdust/300.0)**beta)*DEXP(-1*gamma*(Tdust + T0)/((Tdust**2) + (T0**2)))
                akbar= DEXP(-r(i)%gamma/Tdust)
                IF ( hop_act_competition .EQ. 1 ) THEN
                  ! Revised competition formula
                  akbar = (tunn+akbar*anu0+akbar*anu1)/(tunn+akbar*anu0+akbar*anu1+Rdiff0*&
                    (4.0d0*pi*agr**2.0d0*sitedens)+Rdiff1*(2.0d0*pi*agr**2.0d0*sitedens))
                ELSE
                  ! Just scale exponential
                  akbar = (tunn + (anu0+anu1)*DEXP(-r(i)%gamma/Tdust))/(anu0+anu1)
                  IF ( akbar .GT. 1.0 ) akbar = 1.0
                ENDIF
              ! Special case for H + H2CS ->
              ELSE IF ((s(r(i)%ir1)%name=='gH' .AND. s(r(i)%ir1)%name=='gH2CS') .OR. &
                (s(r(i)%ir2)%name=='gH' .AND. s(r(i)%ir2)%name=='gH2CS')) THEN
                ! -> H2 + HCS
                IF (s(r(i)%ip1)%name=='gH2' .OR. s(r(i)%ip2)%name=='H2') THEN
                  tunn = 1.34e4
                ! -> CH3S
                ELSE IF (s(r(i)%ip1)%name=='gCH3S' .OR. s(r(i)%ip2)%name=='CH3S') THEN
                  tunn = 1.28e9
                ! -> CH2SH
                ELSE IF (s(r(i)%ip1)%name=='gCH2SH' .OR. s(r(i)%ip2)%name=='CH2SH') THEN
                  tunn = 1.54e4
                END IF
                akbar= DEXP(-r(i)%gamma/Tdust)
                IF ( hop_act_competition .EQ. 1 ) THEN
                  ! Revised competition formula
                  akbar = (tunn+akbar*anu0+akbar*anu1)/(tunn+akbar*anu0+akbar*anu1+Rdiff0*&
                    (4.0d0*pi*agr**2.0d0*sitedens)+Rdiff1*(2.0d0*pi*agr**2.0d0*sitedens))
                ELSE
                  ! Just scale exponential
                  akbar = (tunn + (anu0+anu1)*DEXP(-r(i)%gamma/Tdust))/(anu0+anu1)
                  IF ( akbar .GT. 1.0 ) akbar = 1.0
                ENDIF
              ! Special case for H + CH3SH ->
              ELSE IF ((s(r(i)%ir1)%name=='gH' .AND. s(r(i)%ir1)%name=='gCH3SH') .OR. &
                (s(r(i)%ir2)%name=='gH' .AND. s(r(i)%ir2)%name=='gCH3SH')) THEN
                ! -> H2 + CH2SH
                IF (s(r(i)%ip1)%name=='gCH2SH' .OR. s(r(i)%ip2)%name=='CH2SH') THEN
                  tunn = 1.0e8
                ! -> H2 + CH3S
                ELSE IF (s(r(i)%ip1)%name=='gCH3S' .OR. s(r(i)%ip2)%name=='CH3S') THEN
                  tunn = 1.87e2
                ! -> H2S + CH3
                ELSE IF (s(r(i)%ip1)%name=='gCH3' .OR. s(r(i)%ip2)%name=='CH3') THEN
                  tunn = 4.19e2
                END IF
                akbar= DEXP(-r(i)%gamma/Tdust)
                IF ( hop_act_competition .EQ. 1 ) THEN
                  ! Revised competition formula
                  akbar = (tunn+akbar*anu0+akbar*anu1)/(tunn+akbar*anu0+akbar*anu1+Rdiff0*&
                    (4.0d0*pi*agr**2.0d0*sitedens)+Rdiff1*(2.0d0*pi*agr**2.0d0*sitedens))
                ELSE
                  ! Just scale exponential
                  akbar = (tunn + (anu0+anu1)*DEXP(-r(i)%gamma/Tdust))/(anu0+anu1)
                  IF ( akbar .GT. 1.0 ) akbar = 1.0
                ENDIF
              ENDIF
            ENDIF


            r(i)%rate = r(i)%alpha*akbar*(Rdiff0 + Rdiff1)

            IF (r(i)%r1 == r(i)%r2) r(i)%rate = r(i)%rate/2.d0

            IF ( ISNAN(r(i)%rate) .EQV. .TRUE. ) THEN
              r(i)%rate = 0.0
            ENDIF
        CASE (14) !Bulk two-body reaction

            ! Compute their vibrational frequencies:
            IF ( FIXED_NU .EQ. 1 ) THEN
              anu0 = TRIAL_NU
              anu1 = TRIAL_NU
            ELSE
              anu0 = DSQRT(sitedens*ak_B/aMp*s(r(i)%ir1)%edes/PI/PI/s(r(i)%ir1)%weight)
              anu1 = DSQRT(sitedens*ak_B/aMp*s(r(i)%ir2)%edes/PI/PI/s(r(i)%ir2)%weight)
            ENDIF

            ! Compute diffusion rates of the reactants (thermal hopping and quantum tunneling for the model of HHL 1992):
            Rdiff0 = anu0*DEXP(-ebed_factor*ebed*s(r(i)%ir1)%edes/Tdust)
            Rdiff1 = anu1*DEXP(-ebed_factor*ebed*s(r(i)%ir2)%edes/Tdust)

            !Diffusion barriers with tunneling:
            IF (tunneling==1 .AND. FAST_BULK.EQ.0) THEN
              IF (s(r(i)%ir1)%name=='bH' .OR. s(r(i)%ir1)%name=='bH2') THEN
                amur = s(r(i)%ir1)%weight*aMp
                akbar=DEXP(-4D0*PI/hp*1D-8*DSQRT(2D0*amur*ebed_factor*ebed*s(r(i)%ir1)%edes*ak_B))
                Rdiff0 = anu0*akbar
              ENDIF
              IF (s(r(i)%ir2)%name=='bH' .OR. s(r(i)%ir2)%name=='bH2') THEN
                amur = s(r(i)%ir2)%weight*aMp
                akbar=DEXP(-4D0*PI/hp*1D-8*DSQRT(2D0*amur*ebed_factor*ebed*s(r(i)%ir2)%edes*ak_B))
                Rdiff1 = anu1*akbar
              ENDIF

              ! ED = 1040 K = 520 K / 0.5, with 520 K being the best fit barrier height
              ! a = 0.7 Angstrom
              ! See Minissale et al. 2014
              IF (s(r(i)%ir1)%name=='bO') THEN
                amur = s(r(i)%ir2)%weight*aMp
                akbar=DEXP(-4D0*PI/hp*0.7D-8*DSQRT(2D0*amur*ebed_factor*ebed*1040.0*ak_B))
                Rdiff0 = anu0*akbar
              ELSE IF (s(r(i)%ir2)%name=='bO') THEN
                amur = s(r(i)%ir2)%weight*aMp
                akbar=DEXP(-4D0*PI/hp*0.7D-8*DSQRT(2D0*amur*ebed_factor*ebed*1040.0*ak_B))
                Rdiff1 = anu1*akbar
              ENDIF
            ENDIF


            s(r(i)%ir1)%bdiffrate = Rdiff0
            s(r(i)%ir2)%bdiffrate = Rdiff1

            !Reaction activation barriers:
            akbar= DEXP(-r(i)%gamma/Tdust)

            IF (barrier_tunneling==1) THEN
                IF (s(r(i)%ir1)%name=='bH' .OR. s(r(i)%ir1)%name=='bH2' .OR. s(r(i)%ir2)%name=='bH' .OR. s(r(i)%ir2)%name=='bH2') THEN
                    amur = s(r(i)%ir1)%weight*s(r(i)%ir2)%weight/(s(r(i)%ir1)%weight+s(r(i)%ir2)%weight)*aMp
                    akbar=DEXP(-4d0*PI/hp*1E-8*DSQRT(2d0*amur*r(i)%gamma*ak_B))
                ENDIF
            ENDIF

!            IF (hop_act_competition==1) THEN
                ! Ensure that species isn't a radical
!                reactant_names = (/ s(r(i)%ir1)%name , s(r(i)%ir2)%name /)
!                IF ( ANY( radical_names .EQ. TRIM(s(r(i)%ir1)%name) ) .OR. ANY( radical_names .EQ. TRIM(s(r(i)%ir2)%name) ) ) THEN
!                  PRINT *, TRIM(s(r(i)%ir1)%name)," or ",TRIM(s(r(i)%ir2)%name), " should be a radical"
!                ELSE
!                  IF (r(i)%gamma/=0.d0) akbar = 1.0d0 !(akbar*anu0+akbar*anu1)/(akbar*anu0+akbar*anu1+Rdiff0+Rdiff1)
!                ENDIF
!                IF (r(i)%gamma/=0.d0) akbar = 1.0d0 !!!!!!!!!!!!!(akbar*anu0+akbar*anu1)/&
!                   (akbar*anu0+akbar*anu1+Rdiff0*(2.0d0*pi*agr**2.0d0*sitedens)+Rdiff1*&
!                   (2.0d0*pi*agr**2.0d0*sitedens))
!            ENDIF


            IF (shingledecker_tunn==1) THEN
              ! Special case for H2 + OH -> H2O + H
              IF ((s(r(i)%ir1)%name=='bH2' .AND. s(r(i)%ir2)%name=='bOH') .OR. &
                (s(r(i)%ir2)%name=='bH2' .AND. s(r(i)%ir1)%name=='bOH')) THEN
                alpha = 7.64e10
                beta  = 1.00
                gamma = 1339.9
                T0    = 153.2
                tunn  = alpha*((Tdust/300.0)**beta)*DEXP(-1*gamma*(Tdust + T0)/((Tdust**2) + (T0**2)))
                akbar= DEXP(-r(i)%gamma/Tdust)
                ! Scale the activation energy exponential to be the sum of tunneling + thermal
                akbar = (tunn+akbar*(anu0+anu1))/(anu0+anu1)
                IF (akbar .GT. 1.0 ) akbar = 1.0
              ! Special case for H2 + O3 -> H2O + O2
!              ELSE IF ((s(r(i)%ir1)%name=='bH2' .AND. s(r(i)%ir2)%name=='bO3') .OR. &
!                (s(r(i)%ir2)%name=='bH2' .AND. s(r(i)%ir1)%name=='bO3')) THEN
!                alpha = 7.64e10
!                beta  = 1.00
!                gamma = 1339.9
!                T0    = 153.2
!                tunn  = alpha*((Tdust/300.0)**beta)*DEXP(-1*gamma*(Tdust + T0)/((Tdust**2) + (T0**2)))
!                amur  = s(r(i)%ir1)%weight*s(r(i)%ir2)%weight/(s(r(i)%ir1)%weight+s(r(i)%ir2)%weight)*aMp
!                akbar =DEXP(-4d0*PI/hp*barrier_tunneling_w*DSQRT(2d0*amur*r(i)%gamma*ak_B))
!                akbar = (tunn+akbar*anu0+akbar*anu1)/(tunn+akbar*anu0+akbar*anu1+Rdiff0*&
!                  (2.0d0*pi*agr**2.0d0*sitedens)+Rdiff1*(2.0d0*pi*agr**2.0d0*sitedens))
              ! Special case for H + H2O2 -> H2O + OH
              ELSE IF ((s(r(i)%ir1)%name=='bH2O2' .AND. s(r(i)%ir2)%name=='bH') .OR. &
                (s(r(i)%ir2)%name=='bH2O2' .AND. s(r(i)%ir1)%name=='bH')) THEN
                alpha = 1.51e10
                beta  = 0.86
                gamma = 1750.0
                T0    = 180.0
                tunn  = alpha*((Tdust/300.0)**beta)*DEXP(-1*gamma*(Tdust + T0)/((Tdust**2) + (T0**2)))
                akbar= DEXP(-r(i)%gamma/Tdust)
                ! Scale the activation energy exponential to be the sum of tunneling + thermal
                akbar = (tunn+akbar*(anu0+anu1))/(anu0+anu1)
                IF (akbar .GT. 1.0 ) akbar = 1.0
              ! Special case for H + H2S ->  HS + H2
              ELSE IF ((s(r(i)%ir1)%name=='bH' .AND. s(r(i)%ir2)%name=='bH2S') .OR. &
                (s(r(i)%ir2)%name=='bH' .AND. s(r(i)%ir1)%name=='bH2S')) THEN
                alpha = 2.20111E11
                beta  = 0.475667
                gamma = 1399.51
                T0    = 180
                tunn  = alpha*((Tdust/300.0)**beta)*DEXP(-1*gamma*(Tdust + T0)/((Tdust**2) + (T0**2)))
                akbar= DEXP(-r(i)%gamma/Tdust)
                ! Scale the activation energy exponential to be the sum of tunneling + thermal
                akbar = (tunn+akbar*(anu0+anu1))/(anu0+anu1)
                IF (akbar .GT. 1.0 ) akbar = 1.0
              ! Special case for H + CS -> HCS
              ELSE IF ((s(r(i)%ir1)%name=='bH' .AND. s(r(i)%ir1)%name=='bCS') .OR. &
                (s(r(i)%ir2)%name=='bH' .AND. s(r(i)%ir2)%name=='bCS')) THEN
                alpha = 2.05296E12
                beta  = -8.10946
                gamma = 984.986
                T0    = 50.1934
                tunn  = alpha*((Tdust/300.0)**beta)*DEXP(-1*gamma*(Tdust + T0)/((Tdust**2) + (T0**2)))
                akbar= DEXP(-r(i)%gamma/Tdust)
                ! Scale the activation energy exponential to be the sum of tunneling + thermal
                akbar = (tunn+akbar*(anu0+anu1))/(anu0+anu1)
                IF (akbar .GT. 1.0 ) akbar = 1.0
              ! Special case for H + H2CS ->
              ELSE IF ((s(r(i)%ir1)%name=='bH' .AND. s(r(i)%ir1)%name=='bH2CS') .OR. &
                (s(r(i)%ir2)%name=='bH' .AND. s(r(i)%ir2)%name=='bH2CS')) THEN
                ! -> H2 + HCS
                IF (s(r(i)%ip1)%name=='bH2') THEN
                  tunn = 1.34e4
                ! -> CH3S
                ELSE IF (s(r(i)%ip1)%name=='bCH3S') THEN
                  tunn = 1.28e9
                ! -> CH2SH
                ELSE IF (s(r(i)%ip1)%name=='bCH2SH') THEN
                  tunn = 1.54e4
                END IF
                akbar= DEXP(-r(i)%gamma/Tdust)
                ! Scale the activation energy exponential to be the sum of tunneling + thermal
                akbar = (tunn+akbar*(anu0+anu1))/(anu0+anu1)
                IF (akbar .GT. 1.0 ) akbar = 1.0
              ! Special case for H + CH3SH ->
              ELSE IF ((s(r(i)%ir1)%name=='bH' .AND. s(r(i)%ir1)%name=='bCH3SH') .OR. &
                (s(r(i)%ir2)%name=='bH' .AND. s(r(i)%ir2)%name=='bCH3SH')) THEN
                ! -> H2 + CH2SH
                IF (s(r(i)%ip1)%name=='bCH2SH' ) THEN
                  tunn = 1.0e8
                ! -> H2 + CH3S
                ELSE IF (s(r(i)%ip1)%name=='bCH3S') THEN
                  tunn = 1.87e2
                ! -> H2S + CH3
                ELSE IF (s(r(i)%ip1)%name=='bCH3') THEN
                  tunn = 4.19e2
                END IF
                akbar= DEXP(-r(i)%gamma/Tdust)
                ! Scale the activation energy exponential to be the sum of tunneling + thermal
                akbar = (tunn+akbar*(anu0+anu1))/(anu0+anu1)
                IF (akbar .GT. 1.0 ) akbar = 1.0
              ENDIF
            ENDIF

            ! Calculate rate coefficient
            r1 = r(i)%r1
            r2 = r(i)%r2

            ! When FAST_BULK is switched on, reactive species and radicals react
            ! non-diffusively in the bulk of the ice.
            !
            ! Author: C. N. Shingledecker
            IF ( &
                ( FAST_BULK .EQ. 1 ) .AND. &
                ( ANY(radical_names .EQ. r1) .OR.  ANY(radical_names .EQ. r2) ) &
               ) THEN
              r(i)%rate = r(i)%alpha*akbar*(anu0 + anu1)
            ELSE
              r(i)%rate = r(i)%alpha*akbar*(Rdiff0 + Rdiff1)
            ENDIF

!            IF ( ANY(radical_names .EQ. r1) .OR. ANY(radical_names .EQ. r2) ) THEN
!              PRINT *, r_alpha,akbar
!              PRINT *, rate
!            ENDIF

            IF (r(i)%r1 == r(i)%r2) r(i)%rate = r(i)%rate/2.d0

            IF ( ISNAN(r(i)%rate) .EQV. .TRUE. ) THEN
              r(i)%rate = 0.0
            ENDIF

        CASE (15) !Suprathermal surface or bulk reactions
            ! Compute their vibrational frequencies:
            anu0 = DSQRT(sitedens*ak_B/aMp*s(r(i)%ir1)%edes/PI/PI/s(r(i)%ir1)%weight)
            anu1 = DSQRT(sitedens*ak_B/aMp*s(r(i)%ir2)%edes/PI/PI/s(r(i)%ir2)%weight)
            
            IF ( ISNAN(anu0) ) THEN
              PRINT *, "anu0 = NaN"
              anu0 = 1.0e14 ! For electrons. Set characteristic vibration to 1e14 s-1
            ENDIF
 
            IF ( ISNAN(anu1) ) THEN
              PRINT *, "anu1 = NaN"
              anu1 = 1.0e14
            ENDIF
          

            ! Compute reaction rates - no diffusion - of the reactants (Shingledecker et al. 2018):
            IF ( r(i)%r1(1:1) .EQ. 'g' ) THEN
              Rdiff0 = anu0/(4.0d0*pi*agr**2.0d0*sitedens)
              Rdiff1 = anu1/(4.0d0*pi*agr**2.0d0*sitedens)
            ENDIF

            r(i)%rate = r(i)%alpha*(Rdiff0 + Rdiff1)

            IF (r(i)%r1 == r(i)%r2) r(i)%rate = r(i)%rate/2.d0
        CASE (16) !Suprathermal surface or bulk reactions
            ! Compute their vibrational frequencies:
            anu0 = DSQRT(sitedens*ak_B/aMp*s(r(i)%ir1)%edes/PI/PI/s(r(i)%ir1)%weight)
            anu1 = DSQRT(sitedens*ak_B/aMp*s(r(i)%ir2)%edes/PI/PI/s(r(i)%ir2)%weight)

            IF ( ISNAN(anu0) ) THEN
              PRINT *, "anu0 = NaN"
              anu0 = 1.0e14
            ENDIF
 
            IF ( ISNAN(anu1) ) THEN
              PRINT *, "anu1 = NaN"
              anu1 = 1.0e14
            ENDIF


            ! Compute reaction rates - no diffusion - of the reactants (Shingledecker et al. 2018):
            Rdiff0 = anu0
            Rdiff1 = anu1

            r(i)%rate = r(i)%alpha*(Rdiff0 + Rdiff1)

            IF (r(i)%r1 == r(i)%r2) r(i)%rate = r(i)%rate/2.d0
!            PRINT *, r(i)%r1," + ",r(i)%r2," -> ",r(i)%p1," + ",r(i)%p2," + ",r(i)%p3
!            PRINT *,"rate= ",r(i)%rate
!            PRINT *,"********************************************************"
        CASE (17) !Radiolysis of grain species
          IF ( MODEL_EXPERIMENT .EQ. 0 ) THEN
            r(i)%rate = suprathermal*radiolysis*r(i)%alpha*(r(i)%gamma/1.0d2)*8.6*(1.287e-15)*(ZetaCR/1.3d-17)
          ELSE
            r(i)%rate = suprathermal*radiolysis*r(i)%alpha*(r(i)%gamma/1.0d2)*PHI_EXP*Se_EXP
!            r(i)%rate = suprathermal*radiolysis*r(i)%alpha*(r(i)%gamma/1.0d2)*PHI_EXP*Se_EXP*(1.0)*RHO_ICE*(1.0/0.7071067)
            IF ( r(i)%rate .GT. 0.0d0 ) THEN
              PRINT *, r(i)%r1," + IONRAD -> ",r(i)%p1," + ",r(i)%p2," + ",r(i)%p3
              PRINT *, "k_rad=",r(i)%rate
              PRINT *, "************************"
            ENDIF
          ENDIF
        CASE (18) ! Quenching of suprathermal species
          ! Compute their vibrational (trial) frequencies:
          anu0 = DSQRT(sitedens*ak_B/aMp*s(r(i)%ir1)%edes/PI/PI/s(r(i)%ir1)%weight)

            IF ( ISNAN(anu0) ) THEN
              PRINT *, "anu0 = NaN"
              anu0 = 1.0e14
            ENDIF
 
          r(i)%rate = r(i)%alpha*anu0
        CASE(19) ! Photoionization
          ! alpha -> branching fractiona
          ! beta  -> σ, the photoionization cross section
          ! gamma -> δ, the fitting value
          ! PHI_EXP -> ϕ, the photon flux
          ! EXTFAC -> accounts for extinction of photons in the bulk
          ! k = fbr*σ*ϕ*δ
          IF ( FIXED_DVAL .EQ. 1 ) THEN
            IF (s(r(i)%ir1)%name(1:1) =='b') THEN
              r(i)%rate = EXTFAC*PHOTOION*r(i)%alpha*r(i)%beta*PHI_EXP*DVAL
            ELSE
              r(i)%rate = PHOTOION*r(i)%alpha*r(i)%beta*PHI_EXP*DVAL
            ENDIF
          ElSE
            IF (s(r(i)%ir1)%name(1:1) =='b') THEN
              r(i)%rate = EXTFAC*PHOTOION*r(i)%alpha*r(i)%beta*PHI_EXP*r(i)%gamma
            ELSE
              r(i)%rate = PHOTOION*r(i)%alpha*r(i)%beta*PHI_EXP*r(i)%gamma
            ENDIF
          ENDIF
        CASE(20) ! Photoexcitation
          ! alpha -> branching fractiona
          ! beta  -> σ, the photoexcitation cross section
          ! gamma -> δ, the fitting value
          ! PHI_EXP -> ϕ, the photon flux
          ! k = fbr*σ*ϕ*δ
          IF ( FIXED_DVAL .EQ. 1 ) THEN
            IF (s(r(i)%ir1)%name(1:1) =='b') THEN
              r(i)%rate = EXTFAC*PHOTOEXC*r(i)%alpha*r(i)%beta*PHI_EXP*DVAL
            ELSE
              r(i)%rate = PHOTOEXC*r(i)%alpha*r(i)%beta*PHI_EXP*DVAL
            ENDIF
          ElSE
            IF (s(r(i)%ir1)%name(1:1) =='b') THEN
              r(i)%rate = EXTFAC*PHOTOEXC*r(i)%alpha*r(i)%beta*PHI_EXP*r(i)%gamma
            ELSE
              r(i)%rate = PHOTOEXC*r(i)%alpha*r(i)%beta*PHI_EXP*r(i)%gamma
            ENDIF
          ENDIF
        CASE DEFAULT
        END SELECT

ENDDO
!write(22,*)'============='
!stop
!pause
1000 FORMAT(1X,I4,1X,I2,1X,2(A10),10X,5(A10),1pE14.8,1X,F5.2,1X,F8.1)

END SUBROUTINE calc_rates

REAL(wp) FUNCTION get_parameter(time, n_par_steps, time_array, parameter_array)
IMPLICIT NONE
REAL(wp) :: time
REAL(wp), DIMENSION(*) :: time_array, parameter_array
INTEGER  :: n_par_steps, i

get_parameter = 0.d0

IF (time<=time_array(1)) THEN
    get_parameter = parameter_array(1)
    RETURN
ENDIF

IF (time>=time_array(n_par_steps)) THEN
    get_parameter = parameter_array(n_par_steps)
    RETURN
ENDIF

i = 1
DO WHILE (time_array(i)<time)
  i = i + 1
ENDDO

get_parameter = parameter_array(i-1) + (parameter_array(i) - parameter_array(i-1)) / (time_array(i) - time_array(i-1)) * (time - time_array(i-1))

END FUNCTION get_parameter

!------------------------------------------------------------------------------
! Subroutine to compute self- and mutual-shielding of CO and H2:
!------------------------------------------------------------------------------
!SUBROUTINE shield(t,nh2,nco,fh2,fco)
SUBROUTINE shield(t,av,fh2,fco)
IMPLICIT NONE
! Global variables:
REAL(wp) t, nh2, nco, fh2, fco, Av
! Local variables:
INTEGER i, ih2, ico
REAL(wp) b5, x, tth1, tth2, tth3
REAL(wp) th1(52), th2(43), th3(43), tnco(52), tnh2(43), tav(43)
! Common block:
REAL(wp) sigma, delta

sigma = 1.59D+21  !mag.^(-1)/cm^2, NH to Av conversion factor for 0.12 mkm grains
delta = 6.0d-5
!Av = Nh2*2.0D0/sigma
nh2 = av*sigma/2.d0
nco = nh2*delta
! Compute H2 self-shielding (Draine & Bertoldi 1996):
b5 = 3.0D0*dsqrt(t/100.0D0)
x = nh2/5.0D+14

fh2 = 0.965/(1.0D0+x/b5)**2.0D0 + 0.035/DSQRT(1.0D0+x)*DEXP(-8.5D-4*DSQRT(1.0D0+x))
! Compute CO self-shielding (Lee et al. 1996):
OPEN(1, file="Lee_ea_17.txt")
!    rewind 1
READ(1,*)
READ(1,*)
READ(1,*)

DO i = 1, 43
    READ(1,*) tnco(i), th1(i), tnh2(i), th2(i), tav(i), th3(i)
ENDDO

DO i = 44, 52
    READ(1,*) tnco(i), th1(i)
ENDDO

! Search N(CO) position:
DO i = 1, 51
    IF (tnco(i)>=nco) EXIT
ENDDO
i = i - 1
IF (nco<1.0D-08) i = 1
ico = i
tth1 = th1(i) + (th1(i+1) - th1(i)) / (tnco(i+1) - tnco(i)) * (nco - tnco(i))

! Search N(H2) position:
DO i = 1, 42
    IF (tnh2(i)>=nh2) EXIT
ENDDO
i = i - 1
IF (nh2<1.0D-08) i = 1
ih2 = i
tth2 = th2(i) + (th2(i+1) - th2(i)) / (tnh2(i+1) - tnh2(i)) * (nh2 - tnh2(i))

! Search Av position:
DO i = 1, 42
    IF (tav(i)>=Av) EXIT
ENDDO
i=i-1
IF (Av<1.0D-08) i = 1
tth3 = th3(i) + (th3(i+1) - th3(i)) / (tav(i+1) - tav(i)) * (Av - tav(i))

! CO self-shielding:
fco = tth1*tth2*tth3
IF (fco<=0.0D0) fco = 1.0d-30

! Exit:
CLOSE (1)
RETURN
END SUBROUTINE shield


END
