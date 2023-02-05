/*User_IO

    General Description:
    ====================
        User_IO takes care of all the user input requests
---------------------------------------------------------------------*/

#ifndef USER_IO_H
#define USER_IO_H

void Initial_Model(double& Msolar, double& Lsolar, double& Rsolar, double& Ms, double& Ls, 
                   double& Rs, double& Teff, double& X, double& Y, double& Z, char& all_new);

void Change_Model(double& Msolar, double& Lsolar, double& Rsolar, double& Ms, double& Ls, 
                   double& Rs, double& Teff, double& X, double& Y, double& Z, char& all_new);

void Too_Many_Zones(int i, double Msolar, double Ms, double M_r, double Lsolar, double Ls, double L_r, 
                    double r, double Rs, double Rsolar, double Teff, double X, double Y, double Z, 
                    double P_0, double T_0, double rho_0, double kappa_0, double epsilon_0, char rc_flag);

void Final_Results(int i, double Msolar, double Ms, double M_r, double Lsolar, double Ls, double L_r, 
                   double r, double Rs, double Rsolar, double Teff, double X, double Y, double Z, 
                   double P, double T, double rho, double kappa, double epsilon, 
                   double P_0, double T_0, double rho_0, double kappa_0, double epsilon_0, char rc_flag);

#endif
