/*Physics
 
    General Description:
    ====================
        This header file contains the interface function definitions 
        for the physics routines required for the construction of stellar models.
  
---------------------------------------------------------------------*/

#ifndef PHYSICS_H
#define PHYSICS_H

using namespace std;

double PTgradient(double Pm, double P, double Tm, double T);

double Specific_Heat_Ratio();

double Density(double T, double P, double mu, int step_size_condition);

double Opacity(double T, double rho, double X, double Z);

double Optical_Depth_Change(double kappa, double kappam, double rho, double rhom, double r, double rm);

double Nuclear(double T, double rho, double X, double Z);

#endif
