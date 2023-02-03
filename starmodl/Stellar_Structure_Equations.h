/*Stellar_Structure_Equations

    General Description:
    ====================
 
        This header file contains the function interface definitions for the basic equations of 
        stellar structure.  The file also contains a driver function interface definition that selects 
        among the required equations for the Runge Kutta routine.
---------------------------------------------------------------------*/

#ifndef STRUCTUREEQNS_H
#define STRUCTUREEQNS_H

double Structure_Eqns(int i, double r, double X, double Z, double Pm, double Tm, 
                      int step_size_condition, char& rc_flag, double S[], bool& ok);

double dPdr(double M_r, double rho, double r);

double dMdr(double r, double rho);

double dLdr(double r, double rho, double epsilon);

double dTdr(double kappa, double rho, double T, double L_r, double r, double mu, 
            double M_r, double gamma, double dlnPdlnT, char& rc_flag);

#endif


