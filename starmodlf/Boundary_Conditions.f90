MODULE Boundary_Conditions
!
!   General Description:
!   ====================
!       This module contains:
!           (a) the starting conditions for the surface inward integration,
!               assuming that the surface pressure, temperature, and density are all zero.
!           (b) extrapolation formuale to estimate central conditions from the last computed zone
!
!---------------------------------------------------------------------

    USE Constants, ONLY :   dp, pi, two_pi, four_pi_o3, G, a => a_rad, a_o3 => a_rad_o3, c, k => k_B, m_H
    PRIVATE             ::  dp, pi, two_pi, four_pi_o3, G, a, a_o3, c, k, m_H

    CONTAINS
    SUBROUTINE Surface(i, Ms, Ls, rm, X, Z, dr, r, P, T, M_r, L_r, rho, kappa, epsilon, good_surface)

!       General Description:
!       ====================
!           Estimate the temperature and pressure of the outermost zone from the zero boundary condition.
!           Electron scattering and H- ion contributions to the opacity are neglected for simplification.

        USE Composition
        USE Physics
        USE Stellar_Structure_Equations
        USE Zone_Quantities, ONLY   :   Mm, Lm, Pm, Tm, dlnPdlnT, gamma, rc_flag

        IMPLICIT NONE
        INTEGER,    INTENT(IN)      ::  i                                   !zone number
        REAL(dp),   INTENT(IN)      ::  Ms, Ls, rm, X, Z
        REAL(dp),   INTENT(INOUT)   ::  dr
        REAL(dp),   INTENT(OUT)     ::  r, P, T, M_r, L_r, rho, kappa, epsilon
        LOGICAL,    INTENT(OUT)     ::  good_surface

        REAL(dp)                    ::  Y, mu, Aop, tog_bf, XCNO
        REAL(dp),   PARAMETER       ::  g_ff = 1                                    !the free-free Gaunt factor is on the order of unity
        REAL(dp),   PARAMETER       ::  A_bf = 4.34E21, A_ff = 3.68E18              !Bound-free and free-free coefficients
        REAL(dp)                    ::  gamma_ratio, kPadiabatic

        REAL(dp),   PARAMETER       ::  maximum   = 1.0E-8                          !Maximum change in Ms and Ls over surface zone

        INTEGER                     ::  j
        INTEGER,    PARAMETER       ::  j_max = 50

        r = rm + dr

        Y  = Helium(X, Z)
        mu = Mean_Molecular_Weight(X, Y, Z)
        gamma = Specific_Heat_Ratio()
        gamma_ratio = gamma/(gamma - 1)

        j = 0
        outerzone: DO
!           Compute the temperature and pressure for the radiative boundary condition
            rc_flag = "r"
            T = G*Ms*(mu*m_H/(4.25*k))*(1/r - 1/rm)                                 !Eq. (L.2); radiative assumption

            IF (i < 2) THEN
                tog_bf = 0.01                                                       !Assume small value for surface
            ELSE
                tog_bf = 2.82*(rho*(1 + X))**0.2                                    !Taken from Novotny (1973), p. 469
            END IF
            Aop = (A_bf/tog_bf)*Z*(1+X) + A_ff*g_ff*(1-Z)*(1+X)                     !From Eq. (9.22) and (9.23)
            P = SQRT((1/4.25)*(16*pi/3)*(G*Ms/Ls)*(a*c*k/(Aop*mu*m_H)))*T**4.25     !Eq. (L.1)

!           If the zone is convective, recompute the adiabatic temperature and pressure
            dlnPdlnT = PTgradient(Pm, P, Tm, T)
            IF (dlnPdlnT < gamma_ratio .AND. i > 2) THEN
                rc_flag = "c"
                kPadiabatic = Pm/Tm**gamma_ratio
                T = G*Ms*(mu*m_H/(k*gamma_ratio))*(1/r - 1/rm)                      !Eq. (L.3)
                P = kPadiabatic*T**gamma_ratio                                      !Eq. (10.83)
            END IF

!           Compute remaining surface quantities
            rho = Density(T, P, mu)
            IF (rho < 0) THEN
                good_surface = .FALSE.
                EXIT outerzone
            END IF
            kappa   = Opacity(T, rho, X, Z)
            XCNO    = CNO(Z)
            epsilon = Nuclear(T, rho, X, Z)

!           Test to be sure that variations in M_r and L_r are not too large
            M_r = Mm + dMdr(r, rho)*dr
            L_r = Lm + dLdr(r, rho, epsilon)*dr
            IF (ABS((Ms - M_r)/Ms) < maximum .AND. ABS((Ls - L_r)/Ls) < maximum) THEN
                good_surface = .TRUE.
                EXIT outerzone
            END IF

!           If changes in M_r and L_r were too large, repeat with one-half the step size
            j = j + 1
            IF (j > j_max) THEN
                WRITE (*,*) "Unable to converge in SUBROUTINE Surface --- Exiting"
                good_surface = .FALSE.
                EXIT outerzone
            END IF
            dr = dr/2
            r = rm + dr
        END DO outerzone

        IF (.NOT. good_surface) THEN
            WRITE (*,'("The last values obtained by SUBROUTINE Surface were: ", /, &
                    "     M_r = ", ES13.6, "   dM_r/Ms = ", ES13.6, /, &
                    "     L_r = ", ES13.6, "   dL_r/Ls = ", ES13.6)')  M_r, (Ms - M_r)/Ms, L_r, (Ls - L_r)/Ls
        END IF
    END SUBROUTINE Surface
!---------------------------------------------------------------------

    SUBROUTINE Core(M_r, L_r, P, T, X, Z, r, P_0, T_0, rho_0, kappa_0, epsilon_0, rc_flag, dlnPdlnT, good_T)

!       General Description:
!       ====================
!           This routine extrapolates from the inner-most zone to obtain estimates of core conditions in the star

        USE Composition
        USE Physics

        IMPLICIT NONE
        REAL(dp),       INTENT(IN)  ::  M_r, L_r, P, T, X, Z, r
        REAL(dp),       INTENT(OUT) ::  P_0, T_0, rho_0, kappa_0, epsilon_0
        REAL(dp),       INTENT(OUT) ::  dlnPdlnT
        LOGICAL,        INTENT(OUT) ::  good_T
        CHARACTER(1),   INTENT(OUT) ::  rc_flag

        REAL(dp)                    ::  Y, mu, dT
        REAL(dp)                    ::  gamma
        REAL(dp),       PARAMETER   ::  converged = 1.0E-8
        INTEGER                     ::  i
        INTEGER,        PARAMETER   ::  i_max = 50

        rho_0     = M_r/(four_pi_o3*r**3)             !Average density of the central ball
        P_0       = P + (two_pi/3)*G*rho_0**2*r**2    !Central pressure, Eq. (L.4)
        epsilon_0 = L_r/M_r                           !Average energy generation rate of the central ball
        
!       Find core temperature by Newton-Raphson method (including radiation pressure)
        Y   = Helium(X, Z)
        mu  = Mean_Molecular_Weight(X, Y, Z)

        IF (rho_0 > 0) THEN
            i = 0
            T_0 = T
            good_T = .TRUE.
            Find_T_0: DO
                i = i + 1
                dT = -f(T_0)/dfdT(T_0)
                IF (ABS(dT/T_0) < converged) EXIT Find_T_0
                T_0 = T_0 + dT
                IF (i > i_max) THEN
                    WRITE (*,*) "Unable to converge on core temperature in SUBROUTINE Core --- Exiting"
                    good_T = .FALSE.
                    EXIT FIND_T_0
                END IF
            END DO Find_T_0
        ELSE
            T_0 = -T
            good_T = .FALSE.
        END IF

        IF (good_T) THEN
            kappa_0  = Opacity(T_0, rho_0, X, Z)
            dlnPdlnT = PTgradient(P, P_0, T, T_0)
            gamma    = Specific_Heat_Ratio()
            IF (dlnPdlnT < (gamma/(gamma - 1))) THEN
                rc_flag = "c"
            ELSE
                rc_flag = "r"
            END IF
        ELSE
            kappa_0  = -99.9
            dlnPdlnT = -99.9
            rc_flag  = "*"
        END IF

        CONTAINS
            REAL(dp) FUNCTION f(T)
                IMPLICIT NONE
                REAL(dp),   INTENT(IN)  ::  T

                f = rho_0*k*T/(mu*m_H) + a_o3*T**4 - P_0    !f = Ideal Gas Law + Radiation Pressure - core P = 0
            END FUNCTION f

            REAL(dp) FUNCTION dfdT(T)
                IMPLICIT NONE
                REAL(dp),   INTENT(IN)  ::  T

                dfdT = rho_0*k/(mu*m_H) + 4*a_o3*T**3       !df/dT
            END FUNCTION dfdT
    END SUBROUTINE Core    
END MODULE Boundary_Conditions