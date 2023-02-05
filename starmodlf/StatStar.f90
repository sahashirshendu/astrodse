PROGRAM StatStar
!
!   General Description:
!   ====================
!       StatStar computes a ZAMS model using a number of simplifying assumptions about the physics.  The code
!       is not designed to produce precise, research-quality models of ZAMS stars; rather, it is meant to 
!       illustrate many of the fundamental physical ideas discussed in
!
!           "An Introduction to Modern Astrophysics"
!           Bradley W. Carroll and Dale A. Ostlie
!           Second Edition, Addison Wesley,   Copyright 2007.
!
!       StatStar performs an inward integration of the stellar structure equations, given values for the
!       star's mass, luminosity, effective temperature, and composition.
!
!       The simplifying assumptions made here include:
!           (a) A constant composition throughout (characteristic of a ZAMS star).
!           (b) The use of the ideal gas law throughout.
!           (c) The gas is assumed to be completely ionized throughout.
!           (d) Radiation pressure is incorporated.
!           (e) Convection, when present, is taken to be purely adiabatic.
!           (f) Opacity is computed by using approximation formulae:
!               1.  Bound-free processes via Kramer's formula.
!               2.  Free-free processes via Kramer's formula.
!               3.  Electron scattering via Thomson formula.
!               4.  H- ion via fitting function.
!           (g) Surface is assumed to have P = 0, T = 0, rho = 0.
!           (h) Outermost (optically thin) zone is assumed to be radiative.
!           (i) No attempt is made to satisfy the Eddington approximation by
!               adjusting the outside boundary condition.
!
!---------------------------------------------------------------------

    USE Constants, ONLY :   dp
    USE User_IO
    USE Boundary_Conditions
    USE Composition
    USE Physics
    USE Zone_Quantities
    USE ODE_Integrator
    USE Stellar_Structure_Equations

    IMPLICIT NONE
    REAL(dp)                                ::  Msolar, Lsolar, Rsolar, Ms, Ls, Rs, Teff, X, Y, Z, mu
    REAL(dp)                                ::  r, P, T, M_r, L_r
    REAL(dp)                                ::  P_0, T_0, rho_0, epsilon_0, kappa_0
    REAL(dp)                                ::  dr
    REAL(dp),   DIMENSION(n)                ::  PMLT0, PMLT
    REAL(dp),                   PARAMETER   ::  dr_over_r = 1.0E-03         !Initial fractional step size
    
    REAL(dp),                   PARAMETER   ::  M_fraction_limit = 0.01     !Mass fraction stop condition
    REAL(dp),                   PARAMETER   ::  L_fraction_limit = 0.10     !Luminosity stop condition
    REAL(dp),                   PARAMETER   ::  r_fraction_limit = 0.02     !radius stop condition
    INTEGER,                    PARAMETER   ::  maximum_zones = 10000       !Maximum number of zones allowed

    INTEGER,                    PARAMETER   ::  n_surface = 1               !Number of surface boundary zones
    INTEGER                                 ::  i                           !Zone counter
    INTEGER                                 ::  ios                         !I/O status flag
    CHARACTER(1)                            ::  all_new = "Y"               !Select new model parameter

    LOGICAL                                 ::  ok_surface, ok_core, ok_Runge
    LOGICAL,                    PARAMETER   ::  adjust_step_size = .TRUE.   !Allow variable step size

    CHARACTER(26)                           ::  format_table = '(I5, 9ES11.3, 2X, A, F5.1)'
    
!---------------------------------------------------------------------

    New_model: DO
        i = 0       !Initialize zone number

        CALL Initial_Model(Msolar, Lsolar, Rsolar, Ms, Ls, Rs, Teff, X, Y, Z, all_new)
        IF (all_new == "E" .OR. all_new == "e") EXIT New_model

        OPEN(UNIT = 10, FILE = "ZAMSmodel.txt", STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ios)
        IF (ios /= 0) THEN
            WRITE (*,*) "Unable to open output file 'ZAMSmodel.txt' --- terminating calculation"
            EXIT New_model
        END IF

!       Write input data to the output file
        WRITE (10,'(T46,"A ZAMS Stellar Model")')
        WRITE (10,'(T46,"--------------------",/)')
        WRITE (10,'(T46,"M    = ", ES23.16E2, " solar")') Msolar
        WRITE (10,'(T46,"L    = ", ES23.16E2, " solar")') Lsolar
        WRITE (10,'(T46,"R    = ", ES23.16E2, " solar")') Rsolar
        WRITE (10,'(T46,"Teff = ", ES23.16E2, " K")')     Teff
        WRITE (10,'(T46,"X    = ", ES23.16E2)')           X
        WRITE (10,'(T46,"Y    = ", ES23.16E2)')           Y
        WRITE (10,'(T46,"Z    = ", ES23.16E2)')           Z
        WRITE (10,'(//)')

!       Set up previous zone values
        Pm   = 0
        Tm   = 0
        Xm   = X
        Zm   = Z
        rm   = Rs
        taum = 0
        rhom = 0
        kappam = 0
        epsilonm = 0
        dlnPdlnT = 99.9
        rc_flag = "r"

        WRITE (10,'(" zone      r         tau     1-M_r/Ms      L_r         T          P         rho        &
                    &kap        eps    dlnPdlnT")')

        WRITE (10,format_table) i, rm, 0.0, 0.0, Ls, Tm, Pm, rhom, kappam, epsilonm, rc_flag, dlnPdlnT

!       Compute surface zones and step size
        Mm = Ms
        Lm = Ls
        rm = Rs
        dr = -dr_over_r*Rs
        step_size_condition = 0
        Surface_boundary: DO
            i = i + 1

!           Update last zone values
            IF (i > 1) THEN
                Mm = M_r
                Lm = L_r
                rm = r
                Pm = P
                Tm = T
                Xm = X
                Zm = Z
                taum = tau
                rhom = rho
                kappam = kappa
                epsilonm = epsilon
            END IF
        
            CALL Surface(i, Mm, Lm, rm, X, Z, dr, r, P, T, M_r, L_r, rho, kappa, epsilon, ok_surface)
            IF (.NOT. ok_surface) EXIT Surface_boundary

            tau = taum + Optical_Depth_Change(kappa, kappam, rho, rhom, r, rm)
            WRITE (10,format_table) i, r, tau, 1 - M_r/Ms, L_r, T, P, rho, kappa, epsilon, rc_flag, dlnPdlnT

            IF (i == n_surface) EXIT Surface_boundary
        END DO Surface_boundary

        Satisfactory_Surface: IF (ok_surface) THEN
!           Load array of first derivatives to start the general inward integration
            Y        = Helium(X, Z)
            mu       = Mean_Molecular_Weight(X, Y, Z)
            gamma    = Specific_Heat_Ratio()
            dlnPdlnT = PTgradient(Pm, P, Tm, T)
            
            dfdr0(1) = dPdr(M_r, rho, r)
            dfdr0(2) = dMdr(r, rho)
            dfdr0(3) = dLdr(r, rho, epsilon)
            dfdr0(4) = dTdr(kappa, rho, T, L_r, r, mu, M_r, gamma, dlnPdlnT)

!           Main inward integration loop
            Main_Loop: DO
                i = i + 1

!               Update last zone values
                Mm = M_r
                Lm = L_r
                Pm = P
                Tm = T
                Xm = X
                Zm = Z
                rm = r
                taum = tau
                rhom = rho
                kappam = kappa
                epsilonm = epsilon
        
                PMLT0 = (/ Pm, Mm, Lm, Tm /)
                CALL RK_4(n, rm, dr, PMLT0, PMLT, dfdr0, Structure_Eqns, ok_Runge)
                IF (.NOT. ok_Runge) EXIT Main_Loop

!               Results from the current step
                P   = PMLT(1)
                M_r = PMLT(2)
                L_r = PMLT(3)
                T   = PMLT(4)

                tau = taum + Optical_Depth_Change(kappa, kappam, rho, rhom, rm + dr, rm)
                
                WRITE (10,format_table) i, r, tau, 1-M_r/Ms, L_r, T, P, rho, kappa, epsilon, rc_flag, dlnPdlnT
                IF ((M_r/Ms < M_fraction_limit .AND. L_r/Ls < L_fraction_limit .AND. r/Rs < r_fraction_limit) &
                        .OR. T < 0 .OR. P < 0) EXIT Main_Loop
                IF (i > maximum_zones) THEN
                    CALL Too_Many_Zones(i, Msolar, Ms, M_r, Lsolar, Ls, L_r, r, Rs, Rsolar, Teff, X, Y, Z, P_0, T_0, &
                        rho_0, kappa_0, epsilon_0, rc_flag)
                    ok_Runge = .FALSE.
                    EXIT Main_Loop
                END IF
                
!               Is it time to change step size?
                IF (adjust_step_size) THEN
                    SELECT CASE (step_size_condition)
                        CASE(0)
                            IF (M_r < 0.99*Ms) THEN
                                dr = -Rs/100
                                step_size_condition = 1
                            END IF
                        CASE(1)
                            IF (ABS(dr) > 5*r) THEN
                                dr = dr/10
                                step_size_condition = 2
                            END IF
                    END SELECT
                END IF
                r = r + dr
            END DO Main_Loop

            Core_Extrapolation: IF (ok_Runge) THEN
!               Determine core conditions
                i = i + 1
                CALL Core(M_r, L_r, P, T, X, Z, r, P_0, T_0, rho_0, kappa_0, epsilon_0, rc_flag, dlnPdlnT, ok_core)
                IF (.NOT. ok_core) THEN
                    WRITE (*,'(/,"WARNING:  There was a problem with the core extrapolation routine",/)')
                END IF

                tau = tau + Optical_Depth_Change(kappa_0, kappa, rho_0, rho, 0.0_dp, r)
                WRITE (10,format_table) i, 0, tau, 1-M_r/Ms, L_r, T_0, P_0, rho_0, kappa_0, epsilon_0, rc_flag, dlnPdlnT

!               Write initial and final conditions to the screen
                CALL Final_Results(i, Msolar, Ms, M_r, Lsolar, Ls, L_r, r, Rs, Rsolar, Teff, X, Y, Z, &
                        P, T, rho, kappa, epsilon, P_0, T_0, rho_0, kappa_0, epsilon_0, rc_flag)
            END IF Core_Extrapolation
        END IF Satisfactory_Surface

!       Does the user want to compute a new model?
        all_new = "Y"
        CALL Change_Model(Msolar, Lsolar, Rsolar, Ms, Ls, Rs, Teff, X, Y, Z, all_new)
        CLOSE (UNIT = 10, IOSTAT = ios)
        IF (ios /= 0) THEN
            WRITE (*,*) "Unable to close the output file - the new model calculation is being aborted"
            EXIT New_model
        END IF
        IF (all_new == "E") Exit New_model
        WRITE (*, '(//)')
    END DO New_model
END PROGRAM StatStar
