/*Stellar_Structure_Equations

    General Description:
    ====================
 
        This file contains the basic equations of stellar structure.  The file also
        contains a driver function that selects among the required equations for the 
        Runge Kutta routine.
---------------------------------------------------------------------*/

#include <iostream>
#include <cmath>

#include "Constants.h"
#include "Composition.h"
#include "Physics.h"
#include "Stellar_Structure_Equations.h"

double Structure_Eqns(int i, double r, double X, double Z, double Pm, double Tm, 
                      int step_size_condition, char& rc_flag, double S[], bool& ok)
{
    ok = true;

    double P   = S[0];
    double M_r = S[1];
    double L_r = S[2];
    double T   = S[3];

    double Y   = Helium(X, Z);
    double mu  = Mean_Molecular_Weight(X, Y, Z);
    double rho = Density(T, P, mu, step_size_condition);

    if (rho < 0) {
        cout << "Density calculation error in FUNCTION Structure_Eqns" << endl;
        ok = false;}

    double dfdr;
    switch (i) {
        case 0:
            if (ok) dfdr = dPdr(M_r, rho, r);
            else dfdr = 0;
            break;

        case 1:
            if (ok) dfdr = dMdr(r, rho);
            else dfdr = 0;
            break;

        case 2:
            if (ok) {
                double epsilon  = Nuclear(T, rho, X, Z);
                dfdr = dLdr(r, rho, epsilon);}
            else dfdr = 0;
            break;

        case 3:
            if (ok) {
                double kappa    = Opacity(T, rho, X, Z);
                double gamma    = Specific_Heat_Ratio();
                double dlnPdlnT = PTgradient(Pm, P, Tm, T);
                dfdr     = dTdr(kappa, rho, T, L_r, r, mu, M_r, gamma, dlnPdlnT, rc_flag);}
            else dfdr    = 0;
            break;
    }
    return (dfdr);
}

//Hydrostatic Equilibrium
double dPdr(double M_r, double rho, double r)
{
    using ModAstroConstants::G;

    return (-G*M_r*rho/(r*r));                      //Eq. (10.6)
}

//Mass Conservation
double dMdr(double r, double rho)
{
    using ModAstroConstants::four_pi;

    return (four_pi*r*r*rho);                       //Eq. (10.7)
}

//Luminosity Gradient
double dLdr(double r, double rho, double epsilon)
{
    using ModAstroConstants::four_pi;

    return (four_pi*r*r*rho*epsilon);               //Eq. (10.36)
}

//Temperature Gradient
double dTdr(double kappa, double rho, double T, double L_r, double r, double mu, 
            double M_r, double gamma, double dlnPdlnT, char& rc_flag)
{
    using ModAstroConstants::four_pi;
    using ModAstroConstants::four_ac_o3;
    using ModAstroConstants::m_H;
    using ModAstroConstants::k_B;
    using ModAstroConstants::G;
    const double k = k_B;
    
    double gamma_ratio = gamma/(gamma - 1);
    double dTdrfun;
    if (dlnPdlnT > gamma_ratio) {                                       //radiation criterion,   Eq. (10.95)
        dTdrfun = -(kappa*rho/pow(T, 3))*(L_r/(four_pi*r*r))/four_ac_o3;//radiation,             Eq. (10.68)
        rc_flag = 'r';}
    else {
        dTdrfun = -(1/gamma_ratio)*(mu*m_H/k)*(G*M_r/(r*r));            //adiabatic convection,  Eq. (10.89)
        rc_flag = 'c';}

    return (dTdrfun);
}