MODULE Zone_Quantities
!
!   General Description:
!   ====================
!       This module holds required data from the previous and current zones
!
!---------------------------------------------------------------------

    USE Constants, ONLY :   dp

    IMPLICIT NONE
    PRIVATE             ::  dp

!   Previous zone data
    REAL(dp)                                ::  Mm, Lm, rm
    REAL(dp)                                ::  Pm, Tm
    REAL(dp)                                ::  Xm, Zm
    REAL(dp)                                ::  rhom, kappam, taum, epsilonm

!   Current zone data
    REAL(dp)                                ::  rho, kappa, tau, epsilon, gamma, dlnPdlnT
    CHARACTER(1)                            ::  rc_flag

!   Number of stellar structure equations
    INTEGER,                    PARAMETER   ::  n = 4
    
!   Current step size flag
    INTEGER                                 ::  step_size_condition

!   The first derivatives from the stellar structure equations to be used by Runge Kutta routines
    REAL(dp),   DIMENSION(n)                ::  dfdr0
END MODULE Zone_Quantities