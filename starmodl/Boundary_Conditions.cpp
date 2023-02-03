/*Boundary_Conditions
!
!   General Description:
!   ====================
!       This module contains:
!           (a) the starting conditions for the surface inward integration,
!               assuming that the surface pressure, temperature, and density are all zero.
!           (b) extrapolation formuale to estimate central conditions from the last computed zone
!
!---------------------------------------------------------------------*/

#include <cmath>
#include <iostream>
#include <iomanip>

#include "Constants.h"
#include "Composition.h"
#include "Physics.h"
#include "Stellar_Structure_Equations.h"

void Surface(int i, double Ms, double Ls, double Mm, double Lm, double rm, double Pm, double Tm, double X, double Z, 
             double& dr, double& r, double& P, double& T, double& M_r, double& L_r, double& rho, double& kappa, 
             double& epsilon, double& dlnPdlnT, double& gamma, char& rc_flag, int step_size_condition, bool& good_surface)
{
//  General Description:
//  ====================
//      Estimate the temperature and pressure of the outermost zone from the zero boundary condition.
//      Electron scattering and H- ion contributions to the opacity are neglected for simplification.

    using ModAstroConstants::G;
    using ModAstroConstants::m_H;
    using ModAstroConstants::pi;
    using ModAstroConstants::a_rad;
    using ModAstroConstants::c;
    using ModAstroConstants::k_B;

    const double a = a_rad;
    const double k = k_B;
    const double g_ff = 1;              //the free-free Gaunt factor is on the order of unity
    const double A_bf = 4.34e21;        //the bound-free coefficient
    const double A_ff = 3.68e18;        //the free-free coefficient
    const double maximum = 1.0e-08;     //Maximum change in Ms and Ls over surface zone
    const int j_max = 50;

    r = rm + dr;
    double Y  = Helium(X, Z);
    double mu = Mean_Molecular_Weight(X, Y, Z);
    gamma = Specific_Heat_Ratio();
    double gamma_ratio = gamma/(gamma - 1);

    int j = 0;
    bool exit_outerzone = false;
    while (!exit_outerzone){
        //Compute the temperature and pressure for the radiative boundary condition
        rc_flag = 'r';
        T = G*Ms*(mu*m_H/(4.25*k))*(1/r - 1/rm);        //Eq. (L.2); radiative assumption

        double tog_bf;
        if (i < 2) tog_bf = 0.01;                       //Assume small value for surface
        else tog_bf = 2.82*pow((rho*(1 + X)), 0.2);     //Taken from Novotny (1973), p. 469

        double Aop = (A_bf/tog_bf)*Z*(1+X) + A_ff*g_ff*(1-Z)*(1+X);     //From Eq. (9.22) and (9.23)
        P = sqrt((1/4.25)*(16*pi/3)*(G*Ms/Ls)*(a*c*k/(Aop*mu*m_H)))*pow(T, 4.25);    //Eq. (L.1)

        //If the zone is convective, recompute the adiabatic temperature and pressure
        dlnPdlnT = PTgradient(Pm, P, Tm, T);
        double kPadiabatic;
        if (dlnPdlnT < gamma_ratio && i > 2) {
            rc_flag = 'c';
            kPadiabatic = Pm/pow(Tm, gamma_ratio);
            T = G*Ms*(mu*m_H/(k*gamma_ratio))*(1/r - 1/rm);             //Eq. (L.3)
            P = kPadiabatic*pow(T, gamma_ratio);                        //Eq. (10.83)
        }
        
        //Compute remaining surface quantities
        rho = Density(T, P, mu, step_size_condition);
        if (rho < 0) {
            good_surface = false;
            exit_outerzone = true;
        }
        if (exit_outerzone) break;

        kappa = Opacity(T, rho, X, Z);
        double XCNO = CNO(Z);
        epsilon = Nuclear(T, rho, X, Z);

        //Test to be sure that variations in M_r and L_r are not too large
        M_r = Mm + dMdr(r, rho)*dr;
        L_r = Lm + dLdr(r, rho, epsilon)*dr;
        if (fabs((Ms - M_r)/Ms) < maximum && fabs((Ls - L_r)/Ls) < maximum) {
            good_surface = true;
            exit_outerzone = true;
        }
        if (exit_outerzone) break;

        //If changes in M_r and L_r were too large, repeat with one-half the step size
        j++;
        if (j > j_max) {
            cout << "Unable to converge in Function Surface --- Exiting" << endl;
            good_surface = false;
            exit_outerzone = true;
        }
        if (exit_outerzone) break;

        dr /= 2;
        r = rm + dr;
    }

    if (!good_surface) {
        cout << scientific << showpoint << right << setprecision(6)
            << "The last values obtained by Function Surface were: \n"
            << "     M_r = " << setw(13) << M_r
            << "   dM_r/Ms = " << setw(13) << (Ms - M_r)/Ms << "\n"
            << "     L_r = " << setw(13) << L_r
            << "   dL_r/Ls = " << setw(13) << (Ls - L_r)/Ls << endl;
    }
    return;
}

//---------------------------------------------------------------------

void Core(double M_r, double L_r, double P, double T, double X, double Z, double r, double& P_0, double& T_0, 
          double& rho_0, double& kappa_0, double& epsilon_0, char& rc_flag, double& dlnPdlnT, bool& good_T)
{
//  General Description:
//  ====================
//      This routine extrapolates from the inner-most zone to obtain estimates of core conditions in the star

    using ModAstroConstants::four_pi_o3;
    using ModAstroConstants::two_pi;
    using ModAstroConstants::G;
    using ModAstroConstants::m_H;
    using ModAstroConstants::a_rad_o3;
    using ModAstroConstants::k_B;

    const double a_o3 = a_rad_o3;
    const double k = k_B;
    const double converged = 1.0e-08;
    const int i_max = 50;

    rho_0     = M_r/(four_pi_o3*pow(r, 3));         //Average density of the central ball
    P_0       = P + (two_pi/3)*G*rho_0*rho_0*r*r;   //Central pressure, Eq. (L.4)
    epsilon_0 = L_r/M_r;                            //Average energy generation rate of the central ball

    //Find core temperature by Newton-Raphson method (including radiation pressure)
    double Y   = Helium(X, Z);
    double mu  = Mean_Molecular_Weight(X, Y, Z);

    if (rho_0 > 0) {
        int i = 0;
        T_0 = T;
        double dT;
        good_T = true;
        do {
            i++;
            double f = rho_0*k*T_0/(mu*m_H) + a_o3*pow(T_0, 4) - P_0;   //Ideal Gas Law + Radiation Pressure - core P = 0
            double dfdT = rho_0*k/(mu*m_H) + 4*a_o3*pow(T_0, 3);        //df/dT
            
            dT = -f/dfdT;
            T_0 += dT;

            if (i > i_max) {
                cout << "Unable to converge on core temperature in Function Core --- Exiting" << endl;
                good_T = false;}
        } while (i <= i_max && fabs(dT/T_0) >= converged);
    }
    else {
        T_0 = -T;
        good_T = false;}

    if (good_T) {
        kappa_0  = Opacity(T_0, rho_0, X, Z);
        dlnPdlnT = PTgradient(P, P_0, T, T_0);
        double gamma    = Specific_Heat_Ratio();
        if (dlnPdlnT < (gamma/(gamma - 1))) rc_flag = 'c';
        else rc_flag = 'r';}
    else {
        kappa_0  = -99.9;
        dlnPdlnT = -99.9;
        rc_flag  = '*';}

    return;
}