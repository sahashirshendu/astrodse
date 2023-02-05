/*Boundary Conditions

   General Description:
   ====================
       This header file contains function interface definitions for:
           (a) the starting conditions for the surface inward integration,
               assuming that the surface pressure, temperature, and density are all zero.
           (b) extrapolation formuale to estimate central conditions from the last computed zone

           ---------------------------------------------------------------------*/

#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

void Surface(int i, double Ms, double Ls, double Mm, double Lm, double rm, double Pm, double Tm, double X, double Z, 
             double& dr, double& r, double& P, double& T, double& M_r, double& L_r, double& rho, double& kappa, 
             double& epsilon, double& dlnPdlnT, double& gamma, char& rc_flag, int step_size_condition, bool& good_surface);

void Core(double M_r, double L_r, double P, double T, double X, double Z, double r, double& P_0, double& T_0, 
          double& rho_0, double& kappa_0, double& epsilon_0, char& rc_flag, double& dlnPdlnT, bool& good_T);

#endif
