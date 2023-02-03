/*
    General Description:
    ====================
 
        This header file contains the most up-to-date physical and
        astronomical constants in SI units.
 
        "An Introduction to Modern Astrophysics", Appendix I
        Bradley W. Carroll and Dale A. Ostlie
        Addison Wesley, 2007
 
        Weber State University
        Ogden, UT
        modastro@weber.edu
 -------------------------------------------------------------------*/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

namespace ModAstroConstants
{
//The smallest non-zero number and the number of significant figures
    const float         tiny_sp     = 3.4e-38F;
    const double        tiny_dp     = 1.7e-308;
    const long double   tiny_qp     = tiny_dp;
    const short int     sig_fig_sp  = 7;
    const short int     sig_fig_dp  = 15;
    const short int     sig_fig_qp  = sig_fig_dp;
    const float         eps_sp      = 1E-6;
    const double        eps_dp      = 1E-15;
    const long double   eps_qp      = eps_dp;

//The largest number for given precision
    const float         biggest_sp  = 3.4e38F;
    const double        biggest_dp  = 1.7e308;
    const long double   biggest_qp  = biggest_dp;
    const short int     biggest_i2  = 32767;
    const long int      biggest_i4  = 2147483647;
    const long long     biggest_i8  = 9223372036854775807;

//Values related to pi and e
    const long double   pi          = 3.14159265358979323846264338327950;
    const long double   two_pi      = 2*pi;
    const long double   four_pi     = 4*pi;
    const long double   four_pi_o3  = four_pi/3;
    const long double   pi_over_2   = pi/2;
    
    const long double   natural_e   = 2.71828182845904523536028747135266;

//Conversions for radians to degrees and degrees to radians
    const long double   degrees_to_radians = pi/180;
    const long double   radians_to_degrees = 180/pi;

//Physical constants
    const float         G           = 6.673e-11F;
    const double        c           = 2.99792458e08;
    const long double   mu_0        = four_pi*1e-07;
    const long double   epsilon_0   = 1/(mu_0*c*c);

    const double        e_C         = 1.602176462e-19;
    const double        eV          = e_C;
    const double        keV         = eV*1.0e3;
    const double        MeV         = eV*1.0e6;
    const double        GeV         = eV*1.0e9;
    
    const double        h           = 6.62606876e-34;
    const double        hbar        = h/two_pi;

    const double        k_B         = 1.3806503e-23;

    const double        sigma       = 2*pow(pi,5)*pow(k_B,4)/(15*c*c*pow(h,3));
    const double        a_rad       = 4*sigma/c;
    const double        a_rad_o3    = a_rad/3;
    const double        four_ac_o3  = 4*a_rad_o3*c;
    
    const double        m_e         = 9.10938188e-31;
    const double        m_p         = 1.67262158e-27;
    const double        m_n         = 1.67492716e-27;
    const double        m_H         = 1.673532499e-27;
    const double        u           = 1.66053873e-27;
    
    const double        N_A         = 6.02214199e23;
    const double        R_gas       = 8.314472;
    
    const double        a_0         = four_pi*epsilon_0*hbar*hbar/(m_e*e_C*e_C);
    
    const double        R_infty     = m_e*pow(e_C,4)/(64*pow(pi,3)*epsilon_0*epsilon_0*pow(hbar,3)*c);
    const double        R_H         = m_p/(m_e + m_p)*R_infty;

//Time constants
    const short int     hr          = 3600;
    const long int      day         = 24*hr;
    const double        J_yr        = 365.25*day;
    const double        yr          = 3.15581450e7;
    const double        T_yr        = 3.155692519e7;
    const double        G_yr        = 3.1556952e7;
    
//Astronomical length constants
    const double        AU          = 1.4959787066e11;
    const double        pc          = 206264.806*AU;
    const double        ly          = c*J_yr;

//Solar constants
    const float         M_Sun       = 1.9891e30F;
    const float         S_Sun       = 1.365e3F;
    const float         L_Sun       = four_pi*AU*AU*S_Sun;
    const float         R_Sun       = 6.95508e8F;
    const float         Te_Sun      = pow(static_cast <double> (L_Sun/(four_pi*R_Sun*R_Sun*sigma)),0.25);
    
//Solar magnitudes
    const float         Mbol_Sun    =   4.74F;
    const float         MU_Sun      =   5.67F;
    const float         MB_Sun      =   5.47F;
    const float         MV_Sun      =   4.82F;
    const float         Mbol_Sun_ap = -26.83F;
    const float         MU_Sun_ap   = -25.91F;
    const float         MB_Sun_ap   = -26.10F;
    const float         MV_Sun_ap   = -26.75F;
    const float         BC_Sun      =  -0.08F;

//Earth constants
    const float         M_Earth     = 5.9736e24F;
    const float         R_Earth     = 6.378136e6F;
    
//Unit Conversions
    const float         cm          = 1e-2;
    const float         gram        = 1e-3;
    const float         erg         = 1e-7;
    const float         dyne        = 1e-5;
    const double        esu         = 3.335640952e-10;
    const double        statvolt    = 2.997924580e2;
    const float         gauss       = 1e-4;
    const float         angstrom    = 1e-10;
    const float         jansky      = 1e-26;
}

#endif