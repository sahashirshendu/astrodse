MODULE Composition
!
!   General Description:
!   ====================
!       This module contains the information about the composition of the gas
!
!---------------------------------------------------------------------

    USE Constants, ONLY :   dp
    PRIVATE             ::  dp

    CONTAINS
    REAL(dp) FUNCTION Mean_Molecular_Weight(X, Y, Z)   RESULT(mu)

!       General Description:
!       ====================
!           Calculate the mean molecular weight of the gas

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  X, Y, Z

        mu = 1/(2*X + 3*Y/4 + Z/2)           !Assume complete ionization, Eq. (10.16)
    END FUNCTION Mean_Molecular_Weight
!---------------------------------------------------------------------

    REAL(dp) FUNCTION Helium(X, Z) RESULT(Y)

!       General Description:
!       ====================
!           Calculate the amount of Helium-4 in the mixture

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  X, Z
    
        Y = 1 - X - Z
    END FUNCTION Helium
!---------------------------------------------------------------------

    REAL(dp) FUNCTION CNO(Z) RESULT(XCNO)

!       General Description:
!       ====================
!           Calculate the mass fraction of C, N, and O in the mixture

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  Z

        XCNO = Z/2
    END FUNCTION CNO
END MODULE Composition