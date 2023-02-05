/*Physics
 
    General Description:
    ====================
        This module contains the physics routines required for the construction of stellar models.

---------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "Constants.h"
#include "Composition.h"

using namespace std;

double PTgradient(double Pm, double P, double Tm, double T)
{
//      General Description:
//      ====================
//          Compute the pressure gradient with respect to temperature to determine whether convection
//          is required. Limit value of dlnPdlnT for output purposes.

    double dlnPdlnT = ((Tm + T)/(Pm + P))*((Pm - P)/(Tm - T));
    if (dlnPdlnT > 99.9) dlnPdlnT = 99.9;

    return (dlnPdlnT);
}

double Specific_Heat_Ratio()
{
//      General Description:
//      ====================
//          Compute the ratio C_P/C_V
    
    const double monatomic = 5/3.;
        
    return (monatomic);                 // return gamma; assume a purely monatomic gas, Eq. (10.80)
}

double Density(double T, double P, double mu, int step_size_condition)
{
//      General Description:
//      ====================
//          Density computes the density of the gas, assuming the ideal gas law and radiation pressure
//          A negative value for the density indicates that an error was detected in the routine

    using ModAstroConstants::a_rad_o3;
    using ModAstroConstants::m_H;
    using ModAstroConstants::k_B;

    const double a_o3 = a_rad_o3;
    const double k = k_B;
    
    double P_gas = P - a_o3*pow(T, 4);                  //Eq. (10.20)
    if (P_gas <= 0 && T > 0) {                          //Do something desperate
        switch (step_size_condition) {
            case 0:
                P_gas = P;
                break;
            case 1:
                P_gas = 0.001*P;
                break;
            case 2:
                P_gas = 0.0001*P;
                break;
        }
    }

    double rho;
    if (T > 0 && P_gas > 0) rho = P_gas*mu*m_H/(k*T);   //Eq. (10.11)
    else rho = -1;

    if (rho < 0) {
        cout << "A negative density was computed!\n"
            << "Sorry but I am not programmed to handle this new physics :-)\n"
            << "Terminating calculation with: \n"
            << scientific << showpoint << setprecision(6) << right
            << "         T     = " << setw(14) << T << "\n"
            << "         P     = " << setw(14) << P << "\n"
            << "         P_gas = " << setw(14) << P_gas << endl;
    }

    return (rho);
}

double Opacity(double T, double rho, double X, double Z)
{
//      General Description:
//       ====================
//          Opacity computes an approximation of the Rosseland Mean Opacity, based on approximation formulae

    const double A_bf = 4.34e21;
    const double A_ff = 3.68e18;
    const double A_es = 0.02;
    const double A_Hm = 7.9e-34;
    const double g_ff = 1;

    double tog_bf = 0.708*pow(rho*(1 + X), 0.2);                    //Taken from Novotny (1973), p. 469

    double kappa_bf = (A_bf/tog_bf)*Z*(1 + X)*rho/pow(T, 3.5);      //Eq. (9.22)
    double kappa_ff = A_ff*g_ff*(1 - Z)*(1 + X)*rho/pow(T, 3.5);    //Eq. (9.23)
    double kappa_es = A_es*(1 + X);                                 //Eq. (9.27)

    double kappa_Hminus;
    if ((T > 3000 && T < 6000) && (rho > 1e-10 && rho < 1e-5) && (Z > 0.001 && Z < 0.03))
        kappa_Hminus = A_Hm*(Z/0.02)*sqrt(rho)*pow(T, 9);           //Eq. (9.28)
    else kappa_Hminus = 0;

    return (kappa_bf + kappa_ff + kappa_es + kappa_Hminus);         //return the opacity
}

double Optical_Depth_Change(double kappa, double kappam, double rho, double rhom, double r, double rm)
{
//      General Description:
//      ====================
//          Compute the change in optical depth across the zone
    
    return (-(kappa*rho + kappam*rhom)*(r - rm)/2);                 //return dtau; Eq. (9.15)
}

double Nuclear(double T, double rho, double X, double Z)
{
//      General Description:
//      ====================
//          Nuclear computes the nuclear energy generation rates for the proton-proton chains, the CNO cycle,
//          and helium burning.

    const double fpp = 1, f3a = 1;                                  //screening factors
    const double onethird = 1/3.;
    const double twothirds = 2*onethird;
    const double fourthirds = 4*onethird;
    const double fivethirds = 5*onethird;
    const double A_pp = 0.241, A_CNO = 8.67e20, A_He = 50.9;        //reaction rate coefficients

    double T6 = T*1.0e-06;
    double T8 = T*1.0E-08;

    //PP chains (see Hansen and Kawaler, Eq. 6.65, 6.73, and 6.74)
    double psipp = 1 + 1.412e8*(1/X - 1)*exp(-49.98*pow(T6, -onethird));
    double Cpp = 1 + 0.0123*pow(T6, onethird) + 0.0109*pow(T6, twothirds) + 0.000938*T6;
    double eps_pp = A_pp*rho*X*X*fpp*psipp*Cpp*pow(T6, -twothirds)*exp(-33.80*pow(T6, -onethird)); //Eq. (10.46)

    //CNO cycle (Kippenhahn and Weigert, Eq. 18.65)
    double XCNO = CNO(Z);
    double CCNO = 1 + 0.0027*pow(T6, onethird) - 0.00778*pow(T6, twothirds) - 0.000149*T6;
    double eps_CNO = A_CNO*rho*X*XCNO*CCNO*pow(T6, -twothirds)*exp(-152.28*pow(T6, -onethird));     //Eq. (10.58)

    //Helium burning (Kippenhahn and Weigert, Eq. 18.67)
    double Y = Helium(X, Z);
    double eps_He = A_He*rho*rho*pow(Y, 3)/pow(T8, 3)*f3a*exp(-44.027/T8);                          //Eq. (10.62)

    //return epsilon; combined energy generation rate
    return (eps_pp + eps_CNO + eps_He);
}