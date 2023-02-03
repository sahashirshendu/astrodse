/*PROGRAM StatStar
 
    General Description:
    ====================
        StatStar computes a ZAMS model using a number of simplifying assumptions about the physics.  The code
        is not designed to produce precise, research-quality models of ZAMS stars; rather, it is meant to 
        illustrate many of the fundamental physical ideas discussed in
 
            "An Introduction to Modern Astrophysics"
            Bradley W. Carroll and Dale A. Ostlie
            Second Edition, Addison Wesley,   Copyright 2007.
 
        StatStar performs an inward integration of the stellar structure equations, given values for the
        star's mass, luminosity, effective temperature, and composition.
 
        The simplifying assumptions made here include:
            (a) A constant composition throughout (characteristic of a ZAMS star).
            (b) The use of the ideal gas law throughout.
            (c) The gas is assumed to be completely ionized throughout.
            (d) Radiation pressure is incorporated.
            (e) Convection, when present, is taken to be purely adiabatic.
            (f) Opacity is computed by using approximation formulae:
                1.  Bound-free processes via Kramer's formula.
                2.  Free-free processes via Kramer's formula.
                3.  Electron scattering via Thomson formula.
                4.  H- ion via fitting function.
            (g) Surface is assumed to have P = 0, T = 0, rho = 0.
            (h) Outermost (optically thin) zone is assumed to be radiative.
            (i) No attempt is made to satisfy the Eddington approximation by
                adjusting the outside boundary condition.
 
---------------------------------------------------------------------*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "User_IO.h"
#include "Boundary_Conditions.h"
#include "Composition.h"
#include "Physics.h"
#include "ODE_Integrator.h"
#include "Stellar_Structure_Equations.h"

using namespace std;

int main()
{
    const double dr_over_r = 0.001;                 //Initial fractional step size
    const double M_fraction_limit = 0.01;           //Mass fraction stop condition
    const double L_fraction_limit = 0.10;           //Luminosity stop condition
    const double r_fraction_limit = 0.02;           //radius stop condition
    const int    maximum_zones = 10000;             //Maximum number of zones allowed
    const int    n_surface = 1;                     //Number of surface boundary zones
    
    bool ok_surface, ok_core, ok_Runge;
    bool adjust_step_size = true;                   //Allow variable step size
    
    const int n = 4;                                //Number of stellar structure equations
    double dfdr0[n];                               //Stellar structure equation evaluations

    double Mm, Lm, rm, Pm, Tm, Xm, Zm;
    double rhom, kappam, taum, epsilonm;
    double rho, kappa, tau, epsilon, gamma, dlnPdlnT;
    char rc_flag;
    int step_size_condition;

//---------------------------------------------------------------------

    char all_new = 'Y';                             //Select new model parameter

    while (!(all_new == 'E' || all_new == 'e')){

        int i = 0;                                  //Initialize zone number

        double Msolar, Lsolar, Rsolar, Ms, Ls, Rs, Teff, X, Y, Z;
        Initial_Model(Msolar, Lsolar, Rsolar, Ms, Ls, Rs, Teff, X, Y, Z, all_new);
        if (all_new == 'E' || all_new == 'e') exit(0);
        
        //Open the output file
        char xpause;
        ofstream outdata;
        outdata.open("ZAMSmodel.txt", ios::out);
        if (outdata.fail()){
            cout << " Unable to open output file 'ZAMSmodel.txt' --- Terminating calculation"  
                << "\n\n"
                << "Enter any character to quit: ";
            cin  >> xpause;
            exit(1);
        }

        //Write input data to the output file
        outdata 
            << "                                             "
            << "A ZAMS Stellar Model\n" 
            << "                                             "
            << "--------------------\n" << endl;
        outdata << fixed << showpoint << right << setprecision(5)
            << "                                             " << "M    = " << setw(11) << Msolar << " solar" << "\n"
            << "                                             " << "L    = " << setw(11) << Lsolar << " solar" << "\n"
            << "                                             " << "R    = " << setw(11) << Rsolar << " solar" << "\n"
            << "                                             " << "Teff = " << setw(11) << Teff   << " K    " << "\n"
            << "                                             " << "X    = " << setw(11) << X      << "\n"
            << "                                             " << "Y    = " << setw(11) << Y      << "\n"
            << "                                             " << "Z    = " << setw(11) << Z      << "\n\n\n" << endl;

        //Set up previous zone values
        Pm   = 0;
        Tm   = 0;
        Xm   = X;
        Zm   = Z;
        rm   = Rs;
        taum = 0;
        rhom = 0;
        kappam = 0;
        epsilonm = 0;
        dlnPdlnT = 99.9;
        rc_flag = 'r';

        outdata << " zone      r         tau     1-M_r/Ms      L_r         T          P         rho        "
            << "kap        eps    dlnPdlnT" << endl;

        outdata << setw(5) << i
            << scientific << showpoint << right << setprecision(3)
            << setw(11) << rm
            << setw(11) << 0
            << setw(11) << 0
            << setw(11) << Ls
            << setw(11) << Tm
            << setw(11) << Pm
            << setw(11) << rhom
            << setw(11) << kappam
            << setw(11) << epsilonm
            << "  " << rc_flag
            << fixed << showpoint << right << setprecision(1)
            << setw(5) << dlnPdlnT << endl;

        //Compute surface zones and step size
        double r, P, T, M_r, L_r;
        double P_0, T_0, rho_0, epsilon_0, kappa_0;

        Mm = Ms;
        Lm = Ls;
        double dr = -dr_over_r*Rs;
        step_size_condition = 0;
        
        do{
            i++;

            //Update last zone values
            if (i > 1){
                Mm = M_r;
                Lm = L_r;
                rm = r;
                Pm = P;
                Tm = T;
                Xm = X;
                Zm = Z;
                taum = tau;
                rhom = rho;
                kappam = kappa;
                epsilonm = epsilon;
            }
        
            Surface(i, Mm, Lm, Mm, Lm, rm, Pm, Tm, X, Z, dr, r, P, T, M_r, L_r, rho, kappa, epsilon, 
                dlnPdlnT, gamma, rc_flag, step_size_condition, ok_surface);
            if (!ok_surface) break;

            tau = taum + Optical_Depth_Change(kappa, kappam, rho, rhom, r, rm);
            
            outdata << setw(5) << i
                << scientific << showpoint << right << setprecision(3)
                << setw(11) << r
                << setw(11) << tau
                << setw(11) << 1 - M_r/Ms
                << setw(11) << L_r
                << setw(11) << T
                << setw(11) << P
                << setw(11) << rho
                << setw(11) << kappa
                << setw(11) << epsilon
                << "  " << rc_flag
                << fixed << showpoint << right << setprecision(1)
                << setw(5) << dlnPdlnT << endl;
        } while (i < n_surface);

        double mu;
        if (ok_surface){
            //Load array of first derivatives to start the general inward integration
            Y        = Helium(X, Z);
            mu       = Mean_Molecular_Weight(X, Y, Z);
            gamma    = Specific_Heat_Ratio();
            dlnPdlnT = PTgradient(Pm, P, Tm, T);
            
            dfdr0[0] = dPdr(M_r, rho, r);
            dfdr0[1] = dMdr(r, rho);
            dfdr0[2] = dLdr(r, rho, epsilon);
            dfdr0[3] = dTdr(kappa, rho, T, L_r, r, mu, M_r, gamma, dlnPdlnT, rc_flag);

            //Main inward integration loop
            bool exit_main = false;
            while (!exit_main){
                i++;
                
                //Update last zone values
                Mm = M_r;
                Lm = L_r;
                Pm = P;
                Tm = T;
                Xm = X;
                Zm = Z;
                rm = r;
                taum = tau;
                rhom = rho;
                kappam = kappa;
                epsilonm = epsilon;
        
                double PMLT0[4] = {Pm, Mm, Lm, Tm};
                double PMLT[4];
                RK_4(n, rm, dr, PMLT0, PMLT, dfdr0, ok_Runge, X, Z, Pm, Tm, step_size_condition, rc_flag);
                if (!ok_Runge) exit_main = true;
                
                if (ok_Runge){
                    //Results from the current step
                    P   = PMLT[0];
                    M_r = PMLT[1];
                    L_r = PMLT[2];
                    T   = PMLT[3];

                    //Compute current step quantities for output and next step
                    Y   = Helium(X, Z);
                    mu  = Mean_Molecular_Weight(X, Y, Z);
                    rho = Density(T, P, mu, step_size_condition);
                    kappa    = Opacity(T, rho, X, Z);
                    dlnPdlnT = PTgradient(Pm, P, Tm, T);
                    epsilon  = Nuclear(T, rho, X, Z);
                    
                    dfdr0[0] = dPdr(M_r, rho, r);
                    dfdr0[1] = dMdr(r, rho);
                    dfdr0[2] = dLdr(r, rho, epsilon);
                    dfdr0[3] = dTdr(kappa, rho, T, L_r, r, mu, M_r, gamma, dlnPdlnT, rc_flag);
                    
                    tau = taum + Optical_Depth_Change(kappa, kappam, rho, rhom, rm + dr, rm);

                    outdata << setw(5) << i
                        << scientific << showpoint << right << setprecision(3)
                        << setw(11) << r
                        << setw(11) << tau
                        << setw(11) << 1 - M_r/Ms
                        << setw(11) << L_r
                        << setw(11) << T
                        << setw(11) << P
                        << setw(11) << rho
                        << setw(11) << kappa
                        << setw(11) << epsilon
                        << "  "     << rc_flag
                        << fixed    << showpoint << right << setprecision(1)
                        << setw(5)  << dlnPdlnT  << endl;

                    if ((M_r/Ms < M_fraction_limit && L_r/Ls < L_fraction_limit && r/Rs < r_fraction_limit)
                        || T < 0 || P < 0) exit_main = true;
                    else if (i > maximum_zones) {
                        Too_Many_Zones(i, Msolar, Ms, M_r, Lsolar, Ls, L_r, r, Rs, Rsolar, Teff, X, Y, Z, P_0, T_0, 
                            rho_0, kappa_0, epsilon_0, rc_flag);
                        ok_Runge = false;
                        exit_main = true;}

                    if (!exit_main){
                        //Is it time to change step size?
                        if (adjust_step_size){
                            switch (step_size_condition){
                                case 0:
                                    if (M_r < 0.99*Ms) {
                                        dr = -Rs/100;
                                        step_size_condition = 1;}
                                    break;
                                case 1:
                                    if (fabs(dr) > 5*r) {
                                        dr /= 10;
                                        step_size_condition = 2;}
                                    break;
                            }
                        }
                    }
                }
                
                r += dr;
            }

            //Core_Extrapolation
            if (ok_Runge) {
                //Determine core conditions
                i++;
                Core(M_r, L_r, P, T, X, Z, r, P_0, T_0, rho_0, kappa_0, epsilon_0, rc_flag, dlnPdlnT, ok_core);
                if (!ok_core) {
                    cout << "\nWARNING:  There was a problem with the core extrapolation routine" << endl;}

                tau += Optical_Depth_Change(kappa_0, kappa, rho_0, rho, 0, r);
                outdata << setw(5) << i
                    << scientific << showpoint << right << setprecision(3)
                    << setw(11) << 0
                    << setw(11) << tau
                    << setw(11) << 1 - M_r/Ms
                    << setw(11) << L_r
                    << setw(11) << T_0
                    << setw(11) << P_0
                    << setw(11) << rho_0
                    << setw(11) << kappa_0
                    << setw(11) << epsilon_0
                    << "  " << rc_flag
                    << fixed << showpoint << right << setprecision(1)
                    << setw(5) << dlnPdlnT << endl;
                
                //Write initial and final conditions to the screen
                Final_Results(i, Msolar, Ms, M_r, Lsolar, Ls, L_r, r, Rs, Rsolar, Teff, X, Y, Z,
                    P, T, rho, kappa, epsilon, P_0, T_0, rho_0, kappa_0, epsilon_0, rc_flag);
            }
        }

        //Does the user want to compute a new model?
        all_new = 'Y';
        Change_Model(Msolar, Lsolar, Rsolar, Ms, Ls, Rs, Teff, X, Y, Z, all_new);
        outdata.close();
        if (outdata.fail()){
            cout << " Unable to close the output file - the new model calculation is being aborted"  
                << "\n\n"
                << "Enter any character to quit: ";
            cin  >> xpause;
            exit(1);
        }
    }
    return 0;
}