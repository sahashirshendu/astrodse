MODULE Constants
!
!   General Description:
!   ====================
!
!       This module contains the most up-to-date physical and
!       astronomical constants in SI units.  This module also identifies 
!       the correct kind parameters for the current machine.
!
!       "An Introduction to Modern Astrophysics", Appendix I
!       Bradley W. Carroll and Dale A. Ostlie
!       Addison Wesley, 2007
!
!       Weber State University
!       Ogden, UT
!       modastro@weber.edu
!-------------------------------------------------------------------

    IMPLICIT NONE

!The system's precision and range
    INTEGER,    PARAMETER   ::  sp          = SELECTED_REAL_KIND(p = 6, r = 37)
    INTEGER,    PARAMETER   ::  dp          = SELECTED_REAL_KIND(p = 15, r = 307)
    INTEGER,    PARAMETER   ::  qp          = SELECTED_REAL_KIND(p = 33, r = 4931)
    
    INTEGER,    PARAMETER   ::  i1          = SELECTED_INT_KIND(2)
    INTEGER,    PARAMETER   ::  i2          = SELECTED_INT_KIND(4)
    INTEGER,    PARAMETER   ::  i4          = SELECTED_INT_KIND(8)
    INTEGER,    PARAMETER   ::  i8          = SELECTED_INT_KIND(15)

!The smallest non-zero number and the number of significant figures
    REAL(sp),   PARAMETER   ::  tiny_sp     = TINY(1.0_sp)
    REAL(dp),   PARAMETER   ::  tiny_dp     = TINY(1.0_dp)
    REAL(qp),   PARAMETER   ::  tiny_qp     = TINY(1.0_qp)
    INTEGER,    PARAMETER   ::  sig_fig_sp  = PRECISION(1.0_sp)
    INTEGER,    PARAMETER   ::  sig_fig_dp  = PRECISION(1.0_dp)
    INTEGER,    PARAMETER   ::  sig_fig_qp  = PRECISION(1.0_qp)
    REAL(sp),   PARAMETER   ::  eps_sp      = 10.0_sp**(-sig_fig_sp)
    REAL(dp),   PARAMETER   ::  eps_dp      = 10.0_dp**(-sig_fig_dp)
    REAL(qp),   PARAMETER   ::  eps_qp      = 10.0_qp**(-sig_fig_qp)

!The largest number for given precision
    REAL(sp),   PARAMETER   ::  biggest_sp  = HUGE(1.0_sp)
    REAL(dp),   PARAMETER   ::  biggest_dp  = HUGE(1.0_dp)
    REAL(qp),   PARAMETER   ::  biggest_qp  = HUGE(1.0_qp)
    INTEGER(i1),PARAMETER   ::  biggest_i1  = HUGE(1_i1)
    INTEGER(i2),PARAMETER   ::  biggest_i2  = HUGE(1_i2)
    INTEGER(i4),PARAMETER   ::  biggest_i4  = HUGE(1_i4)
    INTEGER(i8),PARAMETER   ::  biggest_i8  = HUGE(1_i8)

!Values related to pi and e
    INTEGER,    PARAMETER,  PRIVATE ::  rpi = SELECTED_REAL_KIND(p = 33, r = 2)
    REAL(rpi),  PARAMETER   ::  pi          = 3.14159265358979323846264338327950_rpi
    REAL(rpi),  PARAMETER   ::  two_pi      = 2*pi
    REAL(rpi),  PARAMETER   ::  four_pi     = 4*pi
    REAL(rpi),  PARAMETER   ::  four_pi_o3  = four_pi/3
    REAL(rpi),  PARAMETER   ::  pi_over_2   = pi/2
    
    REAL(rpi),  PARAMETER   ::  natural_e   = 2.71828182845904523536028747135266_rpi

!Conversions for radians to degrees and degrees to radians
    REAL(rpi),  PARAMETER   ::  degrees_to_radians = pi/180
    REAL(rpi),  PARAMETER   ::  radians_to_degrees = 180/pi

!Physical constants
    INTEGER,    PARAMETER,  PRIVATE :: rG   = SELECTED_REAL_KIND(p = 4, r = 11)
    REAL(rG),   PARAMETER   ::  G           = 6.673e-11_rG
    REAL(qp),   PARAMETER   ::  c           = 2.99792458e08_qp
    REAL(qp),   PARAMETER   ::  mu_0        = four_pi*1e-07_qp
    REAL(qp),   PARAMETER   ::  epsilon_0   = 1/(mu_0*c**2)

    INTEGER,    PARAMETER, PRIVATE  :: reC  = SELECTED_REAL_KIND(p = 10, r = 19)
    REAL(reC),  PARAMETER   ::  e_C         = 1.602176462e-19_reC
    REAL(reC),  PARAMETER   ::  eV          = e_C
    REAL(reC),  PARAMETER   ::  keV         = eV*1.0e3_reC
    REAL(reC),  PARAMETER   ::  MeV         = eV*1.0e6_reC
    REAL(reC),  PARAMETER   ::  GeV         = eV*1.0e9_reC
    
    INTEGER,    PARAMETER,  PRIVATE :: rh   = SELECTED_REAL_KIND(p = 9, r = 34)
    REAL(rh),   PARAMETER   ::  h           = 6.62606876e-34_rh
    REAL(rh),   PARAMETER   ::  hbar        = h/two_pi

    INTEGER,    PARAMETER,  PRIVATE :: rkB  = SELECTED_REAL_KIND(p = 8, r = 23)
    REAL(rkB),  PARAMETER   ::  k_B         = 1.3806503e-23_rkB

    INTEGER,    PARAMETER,  PRIVATE :: rsig = SELECTED_REAL_KIND(p = 8, r = 25)    
    REAL(rsig), PARAMETER   ::  sigma       = 2*pi**5*k_B**4/(15*c**2*h**3)
    REAL(rsig), PARAMETER   ::  a_rad       = 4*sigma/c
    REAL(rsig), PARAMETER   ::  a_rad_o3    = a_rad/3
    REAL(rsig), PARAMETER   ::  four_ac_o3  = 4*a_rad_o3*c
    
    INTEGER,    PARAMETER, PRIVATE  :: rme  = SELECTED_REAL_KIND(p = 9, r = 31)
    REAL(rme),  PARAMETER   ::  m_e         = 9.10938188e-31_rme
    
    INTEGER,    PARAMETER,  PRIVATE :: rmp  = SELECTED_REAL_KIND(p = 9, r = 27)    
    REAL(rmp),  PARAMETER   ::  m_p         = 1.67262158e-27_rmp

    INTEGER,    PARAMETER, PRIVATE  :: rmn  = SELECTED_REAL_KIND(p = 9, r = 27)
    REAL(rmn),  PARAMETER   ::  m_n         = 1.67492716e-27_rmn
    
    INTEGER,    PARAMETER, PRIVATE  :: rmH  = SELECTED_REAL_KIND(p = 10, r = 27)
    REAL(rmH),  PARAMETER   ::  m_H         = 1.673532499e-27_rmH
    
    INTEGER,    PARAMETER, PRIVATE  :: ru   = SELECTED_REAL_KIND(p = 9, r = 27)
    REAL(ru),   PARAMETER   ::  u           = 1.66053873e-27_ru
    
    INTEGER,    PARAMETER, PRIVATE  :: rNA  = SELECTED_REAL_KIND(p = 9, r = 23)
    REAL(rNA),  PARAMETER   ::  N_A         = 6.02214199e23_rNA
    
    INTEGER,    PARAMETER, PRIVATE  :: rR   = SELECTED_REAL_KIND(p = 7, r = 1)
    REAL(rR),   PARAMETER   ::  R_gas       = 8.314472_rR
    
    INTEGER,    PARAMETER, PRIVATE  :: ra0  = SELECTED_REAL_KIND(p = 10, r = 11)
    REAL(ra0),  PARAMETER   ::  a_0         = four_pi*epsilon_0*hbar**2/(m_e*e_C**2)
    
    INTEGER,    PARAMETER, PRIVATE  :: rRH  = SELECTED_REAL_KIND(p = 14, r = 7)
    REAL(rRH),  PARAMETER   ::  R_infty     = m_e*e_C**4/(64*pi**3*epsilon_0**2*hbar**3*c)
    REAL(rRH),  PARAMETER   ::  R_H         = m_p/(m_e + m_p)*R_infty

!Time constants
    INTEGER(i2),PARAMETER   ::  hr          = 3600
    INTEGER(i4),PARAMETER   ::  day         = 24*hr
    REAL(qp),   PARAMETER   ::  J_yr        = 365.25_qp*day

    INTEGER,    PARAMETER, PRIVATE  :: ryr  = SELECTED_REAL_KIND(p = 9, r = 7)
    REAL(ryr),  PARAMETER   ::  yr          = 3.15581450e7_ryr
    
    INTEGER,    PARAMETER, PRIVATE  :: rTyr = SELECTED_REAL_KIND(p = 10, r = 7)
    REAL(rTyr), PARAMETER   ::  T_yr        = 3.155692519e7_rTyr
    
    INTEGER,    PARAMETER, PRIVATE  :: rGyr = SELECTED_REAL_KIND(p = 8, r = 7)
    REAL(rGyr), PARAMETER   ::  G_yr        = 3.1556952e7_rGyr
    
!Astronomical length constants
    INTEGER,    PARAMETER, PRIVATE  :: rAU  = SELECTED_REAL_KIND(p = 11, r = 11)
    REAL(rAU),  PARAMETER   ::  AU          = 1.4959787066e11_dp
    
    INTEGER,    PARAMETER, PRIVATE  :: rpc  = SELECTED_REAL_KIND(p = 11, r = 17)
    REAL(rpc),  PARAMETER   ::  pc          = 206264.806_rpc*AU
    
    REAL(qp),   PARAMETER   ::  ly          = c*J_yr

!Solar constants
    INTEGER,    PARAMETER, PRIVATE  :: rMs  = SELECTED_REAL_KIND(p = 5, r = 30)
    REAL(rMs),  PARAMETER   ::  M_Sun       = 1.9891e30_rMs
    
    INTEGER,    PARAMETER, PRIVATE  :: rSs  = SELECTED_REAL_KIND(p = 4, r = 3)
    REAL(rSs),  PARAMETER   ::  S_Sun       = 1.365e3_rSs
    
    INTEGER,    PARAMETER, PRIVATE  :: rLs  = SELECTED_REAL_KIND(p = 4, r = 26)
    REAL(rLs),  PARAMETER   ::  L_Sun       = four_pi*AU**2*S_Sun
    
    INTEGER,    PARAMETER, PRIVATE  :: rRs  = SELECTED_REAL_KIND(p = 6, r = 6)
    REAL(rMs),  PARAMETER   ::  R_Sun       = 6.95508e8_rRs
    
    INTEGER,    PARAMETER, PRIVATE  :: rTs  = SELECTED_REAL_KIND(p = 4, r = 26)
    REAL(rTs),  PARAMETER   ::  Te_Sun      = (L_Sun/(four_pi*R_Sun**2*sigma))**0.25_qp
    
!Solar magnitudes
    REAL(sp),   PARAMETER   ::  Mbol_Sun    =   4.74
    REAL(sp),   PARAMETER   ::  MU_Sun      =   5.67
    REAL(sp),   PARAMETER   ::  MB_Sun      =   5.47
    REAL(sp),   PARAMETER   ::  MV_Sun      =   4.82
    REAL(sp),   PARAMETER   ::  Mbol_Sun_ap = -26.83
    REAL(sp),   PARAMETER   ::  MU_Sun_ap   = -25.91
    REAL(sp),   PARAMETER   ::  MB_Sun_ap   = -26.10
    REAL(sp),   PARAMETER   ::  MV_Sun_ap   = -26.75
    REAL(sp),   PARAMETER   ::  BC_Sun      =  -0.08

!Earth constants
    INTEGER,    PARAMETER, PRIVATE  :: rMea = SELECTED_REAL_KIND(p = 5, r = 24)
    REAL(rMea), PARAMETER   ::  M_Earth     = 5.9736e24_rMea
    
    INTEGER,    PARAMETER, PRIVATE  :: rRea = SELECTED_REAL_KIND(p = 7, r = 6)
    REAL(rRea), PARAMETER   ::  R_Earth     = 6.378136e6_rRea
    
!Unit Conversions
    REAL(sp),   PARAMETER   ::  cm          = 1e-2
    REAL(sp),   PARAMETER   ::  gram        = 1e-3
    REAL(sp),   PARAMETER   ::  erg         = 1e-7
    REAL(sp),   PARAMETER   ::  dyne        = 1e-5
    REAL(dp),   PARAMETER   ::  esu         = 3.335640952e-10
    REAL(dp),   PARAMETER   ::  statvolt    = 2.997924580e2
    REAL(sp),   PARAMETER   ::  gauss       = 1e-4
    REAL(sp),   PARAMETER   ::  angstrom    = 1e-10
    REAL(sp),   PARAMETER   ::  jansky      = 1e-26
END MODULE Constants