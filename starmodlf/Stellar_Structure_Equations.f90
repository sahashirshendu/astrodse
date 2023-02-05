MODULE Stellar_Structure_Equations

!   General Description:
!   ====================
!
!       This module contains the basic equations of stellar structure.  The module also
!       contains a driver function that selects among the required equations for the 
!       Runge Kutta routines.
!---------------------------------------------------------------------

    USE Constants, ONLY         :   dp, G, four_pi, a => a_rad, four_ac_o3, c, m_H, k => k_B
    USE Zone_Quantities, ONLY   :   n
    PRIVATE                     ::  dp, G, four_pi, a, four_ac_o3, c, m_H, k

    CONTAINS

!   Driver for stellar structure equations
    REAL(dp) FUNCTION Structure_Eqns(i, r, S, ok)   RESULT(dfdr)

        USE Composition
        USE Physics
        USE Zone_Quantities,  X => Xm, Z => Zm

        IMPLICIT NONE
        INTEGER,                    INTENT(IN)  ::  i
        REAL(dp),                   INTENT(IN)  ::  r       !independent variable
        REAL(dp),   DIMENSION(n),   INTENT(IN)  ::  S       !dependent variables
        LOGICAL,                    INTENT(OUT) ::  ok      !Returns .TRUE. if the derivative calculation was successful

        REAL(dp)                                ::  P, M_r, L_r, T
        REAL(dp)                                ::  Y, XCNO, mu

        ok = .TRUE.

        P   = S(1)
        M_r = S(2)
        L_r = S(3)
        T   = S(4)

        Y   = Helium(X, Z)
        mu  = Mean_Molecular_Weight(X, Y, Z)
        rho = Density(T, P, mu)
        IF (rho < 0) THEN
            WRITE (*, '("Density calculation error in FUNCTION Structure_Eqns")')
            ok = .FALSE.
        END IF

        SELECT CASE (i)
            CASE (1)
                IF (ok) THEN
                    dfdr = dPdr(M_r, rho, r)
                ELSE
                    dfdr = 0
                END IF
                dfdr0(1) = dfdr                     !Save result for next zone start

            CASE (2)
                IF (ok) THEN
                    dfdr = dMdr(r, rho)
                ELSE
                    dfdr = 0
                END IF
                dfdr0(2) = dfdr                     !Save result for next zone start

            CASE (3)
                IF (ok) THEN
                    epsilon  = Nuclear(T, rho, X, Z)
                    dfdr     = dLdr(r, rho, epsilon)
                ELSE
                    dfdr     = 0
                END IF
                dfdr0(3) = dfdr                     !Save result for next zone start

            CASE (4)
                IF (ok) THEN
                    kappa    = Opacity(T, rho, X, Z)
                    gamma    = Specific_Heat_Ratio()
                    dlnPdlnT = PTgradient(Pm, P, Tm, T)
                    dfdr     = dTdr(kappa, rho, T, L_r, r, mu, M_r, gamma, dlnPdlnT)
                ELSE
                    dfdr     = 0
                END IF
                dfdr0(4) = dfdr                     !Save result for next zone start
        END SELECT
    END FUNCTION Structure_Eqns
!---------------------------------------------------------------------

!   Hydrostatic Equilibrium
    REAL(dp) FUNCTION dPdr(M_r, rho, r)
        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  M_r, rho, r

        dPdr = -G*M_r*rho/r**2              !Eq. (10.6)
    END FUNCTION dPdr
!---------------------------------------------------------------------

!   Mass Conservation
    REAL(dp) FUNCTION dMdr(r, rho)
        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  r, rho

        dMdr = four_pi*r**2*rho             !Eq. (10.7)
    END FUNCTION dMdr
!---------------------------------------------------------------------

!   Luminosity Gradient
    REAL(dp) FUNCTION dLdr(r, rho, epsilon)
        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  r, rho, epsilon

        dLdr = four_pi*r**2*rho*epsilon     !Eq. (10.36)
    END FUNCTION dLdr
!---------------------------------------------------------------------

!   Temperature Gradient
    REAL(dp) FUNCTION dTdr(kappa, rho, T, L_r, r, mu, M_r, gamma, dlnPdlnT)
        
        USE Zone_Quantities, ONLY   : rc_flag
        
        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  kappa, rho, T, L_r, r, mu, M_r, gamma, dlnPdlnT
        REAL(dp)                ::  gamma_ratio

        gamma_ratio = gamma/(gamma - 1)
        IF (dlnPdlnT > gamma_ratio) THEN                                !radiation criterion,   Eq. (10.95)
            dTdr = -(kappa*rho/T**3)*(L_r/(four_pi*r**2))/four_ac_o3    !radiation,             Eq. (10.68)
            rc_flag = "r"
        ELSE
            dTdr = -(1/gamma_ratio)*(mu*m_H/k)*(G*M_r/r**2)             !adiabatic convection,  Eq. (10.89)
            rc_flag = "c"
        END IF
    END FUNCTION dTdr
END MODULE Stellar_Structure_Equations