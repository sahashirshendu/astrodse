MODULE Physics
!
!   General Description:
!   ====================
!       This module contains the physics routines required for the construction of stellar models.
!
!---------------------------------------------------------------------

    USE Constants, ONLY :   dp, a_o3 => a_rad_o3, k => k_B, m_H
    PRIVATE             ::  dp, a_o3, k, m_H

    CONTAINS
    REAL(dp) FUNCTION PTgradient(Pm, P, Tm, T) RESULT(dlnPdlnT)

!       General Description:
!       ====================
!           Compute the pressure gradient with respect to temperature to determine whether convection
!           is required. Limit value of dlnPdlnT for output purposes.

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  Pm, P, Tm, T

        dlnPdlnT = ((Tm + T)/(Pm + P))*((Pm - P)/(Tm - T))
        IF (dlnPdlnT > 99.9) dlnPdlnT = 99.9
    END FUNCTION PTgradient
!---------------------------------------------------------------------

    REAL(dp) FUNCTION Specific_Heat_Ratio() RESULT(Gamma)

!       General Description:
!       ====================
!           Compute the ratio C_P/C_V

        IMPLICIT NONE
        REAL(dp),   PARAMETER   ::  monatomic = 5/3.0_dp
        
        gamma = monatomic                               !Assume a purely monatomic gas, Eq. (10.80)
    END FUNCTION Specific_Heat_Ratio
!---------------------------------------------------------------------

    REAL(dp) FUNCTION Density(T, P, mu) RESULT(rho)

!       General Description:
!       ====================
!           Density computes the density of the gas, assuming the ideal gas law and radiation pressure
!           A negative value for the density indicates that an error was detected in the routine

        USE Zone_Quantities, ONLY : step_size_condition

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  T, P, mu
        REAL(dp)                ::  P_gas

        P_gas = P - a_o3*T**4                           !Eq. (10.20)
        IF (P_gas <= 0 .AND. T > 0) THEN                !Do something desperate
            SELECT CASE (step_size_condition)
                CASE (0) 
                    P_gas = P
                CASE (1) 
                    P_gas = 0.001*P
                CASE (2) 
                    P_gas = 0.0001*P
            END SELECT
        END IF
        
        IF (T > 0 .AND. P_gas > 0) THEN
            rho = P_gas*mu*m_H/(k*T)                    !Eq. (10.11)
        ELSE
            rho = -1
        END IF
        IF (rho < 0) THEN
            WRITE (*,'("A negative density was computed!", /, &
                     & "Sorry but I am not programmed to handle this new physics :-)", /, &
                     & "Terminating calculation with: ", /, &
                     & "         T     = ", ES13.6, /, &
                     & "         P     = ", ES13.6, /, &
                     & "         P_gas = ", ES13.6)') T, P, P_gas
        END IF
    END FUNCTION Density
!---------------------------------------------------------------------

    REAL(dp) FUNCTION Opacity(T, rho, X, Z) RESULT(kappa)

!       General Description:
!       ====================
!           Opacity computes an approximation of the Rosseland Mean Opacity, based on approximation formulae

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  T, rho, X, Z
        REAL(dp)                ::  kappa_bf, kappa_ff, kappa_es, kappa_Hminus
        REAL(dp)                ::  tog_bf
        REAL(dp),   PARAMETER   ::  g_ff = 1                    !the free-free Gaunt factor is on the order of unity
        REAL(dp),   PARAMETER   ::  A_bf = 4.34E21, A_ff = 3.68E18, A_es = 0.02, A_Hm = 7.9E-34

        tog_bf = 0.708*(rho*(1 + X))**0.2                       !Taken from Novotny (1973), p. 469

        kappa_bf = (A_bf/tog_bf)*Z*(1 + X)*rho/T**3.5           !Eq. (9.22)
        kappa_ff = A_ff*g_ff*(1 - Z)*(1 + X)*rho/T**3.5         !Eq. (9.23)
        kappa_es = A_es*(1 + X)                                 !Eq. (9.27)
        
        IF ((T > 3000 .AND. T < 6000) .AND. (rho > 1E-10 .AND. rho < 1E-5) .AND. (Z > 0.001 .AND. Z < 0.03)) THEN
            kappa_Hminus = A_Hm*(Z/0.02)*SQRT(rho)*T**9         !Eq. (9.28)
        ELSE
            kappa_Hminus = 0
        END IF

        kappa = kappa_bf + kappa_ff + kappa_es + kappa_Hminus
    END FUNCTION Opacity
!---------------------------------------------------------------------

    REAL(dp) FUNCTION Optical_Depth_Change(kappa, kappam, rho, rhom, r, rm) RESULT(dtau)

!       General Description:
!       ====================
!           Compute the change in optical depth across the zone

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  kappa, kappam, rho, rhom, r, rm

        dtau = -(kappa*rho + kappam*rhom)*(r - rm)/2            !Eq. (9.15)
    END FUNCTION Optical_Depth_Change
!---------------------------------------------------------------------

    REAL(dp) FUNCTION Nuclear(T, rho, X, Z) RESULT(epsilon)

!       General Description:
!       ====================
!           Nuclear computes the nuclear energy generation rates for the proton-proton chains, the CNO cycle,
!           and helium burning.

        USE Composition

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  T, rho, X, Z
        REAL(dp)                ::  psipp, Cpp, XCNO, CCNO, Y
        REAL(dp),   PARAMETER   ::  fpp = 1, f3a = 1                                 !screening factors
        REAL(dp)                ::  eps_pp, eps_CNO, eps_He
        REAL(dp)                ::  T6, T8
        REAL(dp),   PARAMETER   ::  onethird = 1/3.0_dp, twothirds = 2*onethird
        REAL(dp),   PARAMETER   ::  fourthirds = 4*onethird, fivethirds = 5*onethird
        REAL(dp),   PARAMETER   ::  A_pp = 0.241, A_CNO = 8.67E20, A_He = 50.9      !reaction rate coefficients

        T6 = T*1.0E-06
        T8 = T*1.0E-08

!       PP chains (see Hansen and Kawaler, Eq. 6.65, 6.73, and 6.74)
        psipp = 1 + 1.412E8*(1/X - 1)*EXP(-49.98*T6**(-onethird))
        Cpp = 1 + 0.0123*T6**onethird + 0.0109*T6**twothirds + 0.000938*T6
        eps_pp = A_pp*rho*X*X*fpp*psipp*Cpp*T6**(-twothirds)*EXP(-33.80*T6**(-onethird))    !Eq. (10.46)

!       CNO cycle (Kippenhahn and Weigert, Eq. 18.65)
        XCNO = CNO(Z)
        CCNO = 1 + 0.0027*T6**onethird - 0.00778*T6**twothirds - 0.000149*T6  
        eps_CNO = A_CNO*rho*X*XCNO*CCNO*T6**(-twothirds)*EXP(-152.28*T6**(-onethird))       !Eq. (10.58)

!       Helium burning (Kippenhahn and Weigert, Eq. 18.67)
        Y = Helium(X, Z)
        eps_He = A_He*rho**2*Y**3/T8**3*f3a*EXP(-44.027/T8)                                 !Eq. (10.62)

!       Combined energy generation rate
        epsilon = eps_pp + eps_CNO + eps_He
    END FUNCTION Nuclear
END MODULE Physics