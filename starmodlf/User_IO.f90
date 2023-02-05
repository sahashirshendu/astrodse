MODULE User_IO

!   General Description:
!   ====================
!       This module takes care of all the user input requests
!---------------------------------------------------------------------

    USE Constants, ONLY :   dp, M_Sun, L_Sun, R_Sun, four_pi, sigma
    PRIVATE             ::  dp, M_Sun, L_Sun, R_Sun, four_pi, sigma

    CONTAINS
    SUBROUTINE Initial_Model(Msolar, Lsolar, Rsolar, Ms, Ls, Rs, Teff, X, Y, Z, all_new)

!       General Description:
!       ====================
!           Gather the initial model data

        USE Composition
        
        IMPLICIT NONE
        REAL(dp),       INTENT(OUT)     ::  Msolar, Lsolar, Rsolar
        REAL(dp),       INTENT(OUT)     ::  Ms, Ls, Rs, Teff, X, Y, Z
        CHARACTER,      INTENT(INOUT)   ::  all_new
        CHARACTER(1)                    ::  y_n
        INTEGER                         ::  ios

!       Write introductory information for user
        WRITE (*,'(   T15, "StatStar is designed to build a ZAMS star")')

        WRITE (*,'(//,T15, "Details of the code are described in:")')
        WRITE (*,'(   T15, "   An Introduction to Modern Astrophysics")')
        WRITE (*,'(   T15, "   Bradley W. Carroll and Dale A. Ostlie")')
        WRITE (*,'(   T15, "      Second Edition, Addison Wesley")')
        WRITE (*,'(   T15, "      copyright 2007",//)')

        WRITE (*,'(   T15, "The user will be asked to enter the following quantities:")')
        WRITE (*,'(   T15, "   Mass of the star       (in solar units)")')
        WRITE (*,'(   T15, "   Luminosity of the star (in solar units)")')
        WRITE (*,'(   T15, "   Effective Temperature  (in K)")')
        WRITE (*,'(   T15, "   Hydrogen mass fraction (X)")')
        WRITE (*,'(   T15, "   Metals mass fraction   (Z)",//)')

        All_parameters: DO
            IF (all_new == "Y" .OR. all_new == "y" .OR. all_new == "A" .OR. all_new == "a") THEN
                M_loop: DO
                    WRITE (*,'("Enter the mass (in solar units)       :  Msolar = ")', ADVANCE = 'NO')
                    READ  (*,*, iostat = ios) Msolar
                    IF (ios == 0 .AND. Msolar > 0) EXIT M_loop
                    WRITE (*,'("Invalid value entered - please try again",/)')
                END DO M_loop
            
                T_loop: DO
                    WRITE (*,'("Enter the effective temperature (in K):  Teff   = ")', ADVANCE = 'NO') 
                    READ  (*,*, iostat = ios) Teff
                    IF (ios == 0 .AND. Teff > 0) EXIT T_loop
                    WRITE (*,'("Invalid value entered - please try again",/)')
                END DO T_loop

                L_loop: DO
                    WRITE (*,'("Enter the luminosity (in solar units) :  Lsolar = ")', ADVANCE = 'NO') 
                    READ  (*,*, iostat = ios) Lsolar
                    IF (ios == 0 .AND. Lsolar > 0) EXIT L_loop
                    WRITE (*,'("Invalid value entered - please try again",/)')
                END DO L_loop

                XYZ: DO
                    X_loop: DO
                        WRITE (*,'("Enter the mass fraction of hydrogen   :  X      = ")', ADVANCE = 'NO') 
                        READ  (*,*, iostat = ios) X
                        IF (ios == 0) THEN
                            IF (X >= 0 .AND. X <= 1) THEN
                                EXIT X_loop
                            ELSE
                                WRITE (*,'(/,"0 <= X <= 1 is required")')
                            END IF
                        ELSE
                            WRITE (*,'(/)')
                        END IF
                        WRITE (*,'("Invalid value entered - please try again",/)')
                    END DO X_loop

                    Z_loop: DO
                        WRITE (*,'("Enter the mass fraction of metals     :  Z      = ")', ADVANCE = 'NO')
                        READ  (*,*, iostat = ios) Z
                        IF (ios == 0) THEN
                            IF (Z >= 0 .AND. Z <= 1) THEN
                                EXIT Z_loop
                            ELSE
                                WRITE (*,'(/,"0 <= Z <= 1 is required")')
                            END IF
                        ELSE
                            WRITE (*,'(/)')
                        END IF
                        WRITE (*,'("Invalid value entered - please try again",/)')
                    END DO Z_loop
                
                    Y = Helium(X, Z)
                    IF (Y < 0 .OR. Y > 1) THEN
                        WRITE (*,'("Note that 0 <= X + Z <= 1 is required", /, "  Please reenter composition values",/)')
                    ELSE
                        EXIT XYZ
                    END IF
                END DO XYZ
            END IF

!           Compute SI values
            Ms = Msolar*M_Sun
            Ls = Lsolar*L_Sun
            Rs = SQRT(Ls/(four_pi*sigma*Teff**4))   !Eq. (3.17)
            Rsolar = Rs/R_Sun                       !Solar radius from SI value

!           Allow the user the opportunity to change values as needed
            Fix_parameters: DO
                WRITE (*,'(/,"Your model parameters are:", &
                    &      T42, "M      = ", ES23.16E2, " solar")') Msolar
                WRITE (*,'(T42, "Teff   = ", ES23.16E2, " K    ")') Teff
                WRITE (*,'(T42, "L      = ", ES23.16E2, " solar")') Lsolar
                WRITE (*,'(T42, "R      = ", ES23.16E2, " solar")') Rsolar
                WRITE (*,'(T42, "X      = ", ES23.16E2)') X
                WRITE (*,'(T42, "Y      = ", ES23.16E2)') Y
                WRITE (*,'(T42, "Z      = ", ES23.16E2)') Z
                
                WRITE (*,'(/,"Are these values ok (y/n)? ")', ADVANCE = 'NO')
                Yes_No: DO
                    READ  (*,*) y_n
                    IF (y_n == "Y" .OR. y_n == "y" .OR. y_n == "N" .OR. y_n == "n") EXIT Yes_No
                    WRITE (*,'("Please answer Yes (y) or No (n):  ")', ADVANCE = 'NO')
                END DO Yes_No
                IF (y_n == "Y" .OR. y_n == "y") EXIT All_parameters
                
                all_new = "N"
                CALL Change_Model(Msolar, Lsolar, Rsolar, Ms, Ls, Rs, Teff, X, Y, Z, all_new)
                IF (all_new == "E" .OR. all_new == "e") EXIT All_parameters
                IF (all_new == "A" .OR. all_new == "a") EXIT Fix_parameters
            END DO Fix_parameters
        END DO All_parameters
    END SUBROUTINE Initial_Model
!---------------------------------------------------------------------

    SUBROUTINE Change_Model(Msolar, Lsolar, Rsolar, Ms, Ls, Rs, Teff, X, Y, Z, all_new)

!       General Description:
!       ====================
!           Get updated model input data

        USE Composition
        
        IMPLICIT NONE
        REAL(dp),       INTENT(INOUT)   ::  Msolar, Lsolar, Rsolar
        REAL(dp),       INTENT(INOUT)   ::  Ms, Ls, Rs, Teff, X, Y, Z
        CHARACTER(1),   INTENT(INOUT)   ::  all_new
        CHARACTER(1)                    ::  y_n
        INTEGER                         ::  ios

        IF (all_new == "Y" .OR. all_new == "y") THEN
            WRITE (*,'(/,"Would you like to run another model?", /, &
                    &    "Your previous results will be overwritten in the output file. (y/n): ")', ADVANCE = 'NO')
            Yes_No: DO
                READ  (*,*) y_n
                IF (y_n == "Y" .OR. y_n == "y" .OR. y_n == "N" .OR. y_n == "n") EXIT Yes_No
                WRITE (*,'("Please answer Yes (y) or No (n):  ")', ADVANCE = 'NO')
            END DO Yes_No

            IF (y_n == "Y" .OR. y_n == "y") THEN
                WRITE (*, '(/,"Which variable would you like to change?")')
                WRITE (*, '("     M = Mass",                       T40, "Current value = ", ES23.16E2, " solar")') Msolar
                WRITE (*, '("     T = effective Temperature",      T40, "Current value = ", ES23.16E2, " K    ")') Teff
                WRITE (*, '("     L = Luminosity",                 T40, "Current value = ", ES23.16E2, " solar")') Lsolar
                WRITE (*, '("     X = hydrogen mass mraction (X)", T40, "Current value = ", ES23.16E2, "      ")') X
                WRITE (*, '("     Z = metals mass fraction (Z)",   T40, "Current value = ", ES23.16E2, "      ")') Z
                WRITE (*, '("     A = select an All new set of model parameters")')
                WRITE (*, '("     E = Exit the calculation")')
                WRITE (*, '(/, "Select a letter: ")', ADVANCE = 'NO')
                All_New1: DO
                    READ  (*,*) all_new
                    IF (all_new == "M" .OR. all_new == "m" .OR. &
                        all_new == "T" .OR. all_new == "t" .OR. &
                        all_new == "L" .OR. all_new == "l" .OR. &
                        all_new == "X" .OR. all_new == "x" .OR. &
                        all_new == "Z" .OR. all_new == "z" .OR. &
                        all_new == "A" .OR. all_new == "a" .OR. &
                        all_new == "E" .OR. all_new == "e") EXIT All_New1
                    WRITE (*,'("Please respond with one of the options listed:  ")', ADVANCE = 'NO')
                END DO All_New1
            END IF    
        ELSE
            y_n = "Y"
            WRITE (*, '(/,"Which variable would you like to change?")')
            WRITE (*, '("     M = Mass")')
            WRITE (*, '("     T = effective Temperature")')
            WRITE (*, '("     L = Luminosity")')
            WRITE (*, '("     X = hydrogen mass fraction (X)")')
            WRITE (*, '("     Z = metals mass fraction (Z)")')
            WRITE (*, '("     A = select an All new set of model parameters")')
            WRITE (*, '("     E = Exit the calculation")')
            WRITE (*, '(/, "Select a letter: ")', ADVANCE = 'NO')
            All_New2: DO
                READ  (*,*) all_new
                IF (all_new == "M" .OR. all_new == "m" .OR. &
                    all_new == "T" .OR. all_new == "t" .OR. &
                    all_new == "L" .OR. all_new == "l" .OR. &
                    all_new == "X" .OR. all_new == "x" .OR. &
                    all_new == "Z" .OR. all_new == "z" .OR. &
                    all_new == "A" .OR. all_new == "a" .OR. &
                    all_new == "E" .OR. all_new == "e") EXIT All_New2
                WRITE (*,'("Please respond with one of the options listed:  ")', ADVANCE = 'NO')
            END DO All_New2
            IF (all_new == "E" .OR. all_new == "e") y_n = "N"
        END IF

        IF (y_n == "Y" .OR. y_n == "y") THEN
            SELECT CASE (all_new)
                CASE ("M", "m")
                    M_loop: DO
                        WRITE (*,'("Enter the mass (in solar units)       :  Msolar = ")', ADVANCE = 'NO') 
                        READ  (*,*, iostat = ios) Msolar
                        IF (ios == 0 .AND. Msolar > 0) THEN
                            Ms = Msolar*M_Sun
                            EXIT M_loop
                        END IF
                        WRITE (*,'("Invalid value entered - please try again",/)')
                    END DO M_loop
                    all_new = "n"
                CASE ("T", "t")
                    T_loop: DO
                        WRITE (*,'("Enter the effective temperature (in K):  Teff   = ")', ADVANCE = 'NO') 
                        READ  (*,*, iostat = ios) Teff
                        IF (ios == 0 .AND. Teff > 0) THEN
                            Rs = SQRT(Ls/(four_pi*sigma*Teff**4))   !Eq. (3.17)
                            Rsolar = Rs/R_Sun                       !Solar radius from SI value
                            EXIT T_loop
                        END IF
                        WRITE (*,'("Invalid value entered - please try again",/)')
                    END DO T_loop
                    all_new = "n"
                CASE ("L", "l")
                    L_loop: DO
                        WRITE (*,'("Enter the luminosity (in solar units) :  Lsolar = ")', ADVANCE = 'NO') 
                        READ  (*,*, iostat = ios) Lsolar
                        IF (ios == 0 .AND. Lsolar > 0) THEN
                            Ls = Lsolar*L_Sun
                            Rs = SQRT(Ls/(four_pi*sigma*Teff**4))   !Eq. (3.17)
                            Rsolar = Rs/R_Sun                       !Solar radius from SI value
                            EXIT L_loop
                        END IF
                        WRITE (*,'("Invalid value entered - please try again",/)')
                    END DO L_loop
                    all_new = "n"
                CASE ("X", "x")
                    X_loop: DO
                        WRITE (*,'("Enter the mass fraction of hydrogen   :  X      = ")', ADVANCE = 'NO') 
                        READ  (*,*, iostat = ios) X
                        IF (ios == 0) THEN
                            IF (X >= 0 .AND. X <= 1) THEN
                                EXIT X_loop
                            ELSE
                                WRITE (*,'(/,"0 <= X <= 1 is required")')
                            END IF
                        ELSE
                            WRITE (*,'(/)')
                        END IF
                        WRITE (*,'("Invalid value entered - please try again",/)')
                    END DO X_loop
                    X_again: DO
                        Y = Helium(X, Z)
                        IF (Y >= 0 .AND. Y <= 1) EXIT X_again
                        WRITE (*,'("Note that 0 <= X + Z < = 1 is required")')
                        X_repeat: DO
                            WRITE (*,'("  Please try again:                      X      = ")', ADVANCE = 'NO')
                            READ (*,*, iostat = ios) X
                            IF (ios == 0) THEN
                                IF (X >= 0 .AND. X <= 1) THEN
                                    EXIT X_repeat
                                ELSE
                                    WRITE (*,'(/,"0 <= X <= 1 is required")')
                                END IF
                            ELSE
                                WRITE (*,'(/)')
                            END IF
                            WRITE (*,'("Invalid value entered",/)')
                        END DO X_repeat
                    END DO X_again
                    all_new = "n"
                CASE ("Z", "z")
                    Z_loop: DO
                        WRITE (*,'("Enter the mass fraction of metals     :  Z      = ")', ADVANCE = 'NO') 
                        READ  (*,*, iostat = ios) Z
                        IF (ios == 0) THEN
                            IF (Z >= 0 .AND. Z <= 1) THEN
                                EXIT Z_loop
                            ELSE
                                WRITE (*,'(/,"0 <= Z <= 1 is required")')
                            END IF
                        ELSE
                            WRITE (*,'(/)')
                        END IF
                        WRITE (*,'("Invalid value entered - please try again",/)')
                    END DO Z_loop
                    Z_again: DO
                        Y = Helium(X, Z)
                        IF (Y >= 0 .AND. Y <= 1) EXIT Z_again
                        WRITE (*,'("Note that 0 <= X + Z < = 1 is required")')
                        Z_repeat: DO
                            WRITE (*,'("  Please try again:                      Z      = ")', ADVANCE = 'NO')
                            READ (*,*, iostat = ios) Z
                            IF (ios == 0) THEN
                                IF (Z >= 0 .AND. Z <= 1) THEN
                                    EXIT Z_repeat
                                ELSE
                                    WRITE (*,'(/,"0 <= Z <= 1 is required")')
                                END IF
                            ELSE
                                WRITE (*,'(/)')
                            END IF
                            WRITE (*,'("Invalid value entered",/)')
                        END DO Z_repeat
                    END DO Z_again
                    all_new = "n"
                CASE ("E", "e")
                    y_n = "N"
                CASE ("A", "a")
                    y_n = "Y"
                CASE DEFAULT
                    all_new = "y"
            END SELECT
        ELSE
            all_new = "E"           !Exit calculations
        END IF
    END SUBROUTINE Change_Model
!---------------------------------------------------------------------

    SUBROUTINE Too_Many_Zones(i, Msolar, Ms, M_r, Lsolar, Ls, L_r, r, Rs, Rsolar, Teff, X, Y, Z, P_0, T_0, &
                  & rho_0, kappa_0, epsilon_0, rc_flag)

!       General Description:
!       ====================
!           Tell user that the maximum number of zones has been exceeded

        IMPLICIT NONE
        INTEGER,        INTENT(IN)  ::  i
        REAL(dp),       INTENT(IN)  ::  Msolar, Lsolar, Rsolar
        REAL(dp),       INTENT(IN)  ::  Rs, Ms, Ls, Teff, X, Y, Z
        REAL(dp),       INTENT(IN)  ::  r, M_r, L_r
        REAL(dp),       INTENT(IN)  ::  P_0, T_0, rho_0, kappa_0, epsilon_0
        CHARACTER(1),   INTENT(IN)  ::  rc_flag

        WRITE (*,'(/,"The maximum number of zones has been exceeded for this model - Sorry!")')
        WRITE (*,'("The conditions at the time the model was terminated were: ")')
        WRITE (*,'(T5,"Surface Conditions:",           T35, "Last Zone Calculated:")')
        WRITE (*,'(T5,"-------------------",           T35, "---------------------")')
        WRITE (*,'(T5,"M    = ", F12.6, " solar",      T35, "M_r/Ms  = ", F12.6)')                 Msolar, M_r/Ms
        WRITE (*,'(T5,"Teff = ", F12.6, " K",          T35, "L_r/LS  = ", F12.6)')                 Teff,   L_r/Ls
        WRITE (*,'(T5,"L    = ", F12.6, " solar",      T35, "r/Rs    = ", F12.6)')                 Lsolar, r/Rs
        WRITE (*,'(T5,"R    = ", F12.6, " solar",      T35, "P       = ", ES12.5, " N/m^2")')      Rsolar, P_0
        WRITE (*,'(T5,"X    = ", F12.6,                T35, "T       = ", ES12.5, " K")')          X,      T_0
        WRITE (*,'(T5,"Y    = ", F12.6,                T35, "rho     = ", ES12.5, " kg/m^3")')     Y,      rho_0
        WRITE (*,'(T5,"Z    = ", F12.6,                T35, "kappa   = ", ES12.5, " m^2/kg")')     Z,      kappa_0
        WRITE (*,'(                                    T35, "epsilon = ", ES12.5, " W/kg")')               epsilon_0
        IF (rc_flag == "r") THEN
            WRITE (*,'(T35, "The core is RADIATIVE")')
        ELSE
            WRITE (*,'(T35, "The core is CONVECTIVE")')
        END IF
    END SUBROUTINE Too_Many_Zones
!---------------------------------------------------------------------

    SUBROUTINE Final_Results(i, Msolar, Ms, M_r, Lsolar, Ls, L_r, r, Rs, Rsolar, Teff, X, Y, Z, &
                 & P, T, rho, kappa, epsilon, P_0, T_0, rho_0, kappa_0, epsilon_0, rc_flag)

!       General Description:
!       ====================
!           Tell the user the conditions at the surface and core of the completed model

        IMPLICIT NONE
        INTEGER,        INTENT(IN)  ::  i
        REAL(dp),       INTENT(IN)  ::  Msolar, Lsolar, Rsolar
        REAL(dp),       INTENT(IN)  ::  Ms, Ls, Rs, Teff, X, Y, Z
        REAL(dp),       INTENT(IN)  ::  r, M_r, L_r
        REAL(dp),       INTENT(IN)  ::  p, T, rho, kappa, epsilon
        REAL(dp),       INTENT(IN)  ::  P_0, T_0, rho_0, kappa_0, epsilon_0
        CHARACTER(1),   INTENT(IN)  ::  rc_flag

        WRITE (*,'(//,"*********************THE INTEGRATION HAS BEEN COMPLETED*********************")')
        WRITE (*,'(T5,"Surface Conditions:",           T35, "Core Conditions:")')
        WRITE (*,'(T5,"-------------------",           T35, "----------------")')
        WRITE (*,'(T5,"M    = ", F12.6, " solar",      T35, "M_r/Ms  = ", F12.6)')                 Msolar, M_r/Ms
        WRITE (*,'(T5,"Teff = ", F12.6, " K",          T35, "L_r/LS  = ", F12.6)')                 Teff,   L_r/Ls
        WRITE (*,'(T5,"L    = ", F12.6, " solar",      T35, "r/Rs    = ", F12.6)')                 Lsolar, r/Rs
        WRITE (*,'(T5,"R    = ", F12.6, " solar",      T35, "P       = ", ES12.5, " N/m^2")')      Rsolar, P_0
        WRITE (*,'(T5,"X    = ", F12.6,                T35, "T       = ", ES12.5, " K")')          X,      T_0
        WRITE (*,'(T5,"Y    = ", F12.6,                T35, "rho     = ", ES12.5, " kg/m^3")')     Y,      rho_0
        WRITE (*,'(T5,"Z    = ", F12.6,                T35, "kappa   = ", ES12.5, " m^2/kg")')     Z,      kappa_0
        WRITE (*,'(                                    T35, "epsilon = ", ES12.5, " W/kg")')               epsilon_0
        IF (rc_flag == "r") THEN
            WRITE (*,'(T35, "The core is RADIATIVE")')
        ELSE
            WRITE (*,'(T35, "The core is CONVECTIVE")')
        END IF

        WRITE (*,'(/,"For your information, the conditions in the last zone above the core are:")')
        WRITE (*,'(T35, "P       = ", ES12.5, " N/m^2")') P
        WRITE (*,'(T35, "T       = ", ES12.5, " K")') T
        WRITE (*,'(T35, "rho     = ", ES12.5, " kg/m^3")') rho
        WRITE (*,'(T35, "kappa   = ", ES12.5, " m^2/kg")') kappa
        WRITE (*,'(T35, "epsilon = ", ES12.5, " W/kg")') epsilon

        WRITE (*,'(/,"The number of mass shells in this model: ", I5)') i
        WRITE (*,'("The details of the model are available in ZAMSmodel.txt")')
    END SUBROUTINE Final_Results
END MODULE User_IO
