/*User_IO

    General Description:
    ====================
        This file takes care of all the user input requests and most output
---------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "Constants.h"
#include "Composition.h"
#include "User_IO.h"

using namespace std;

void Initial_Model(double& Msolar, double& Lsolar, double& Rsolar, double& Ms, double& Ls, 
                   double& Rs, double& Teff, double& X, double& Y, double& Z, char& all_new)
{
//      General Description:
//      ====================
//          Gather the initial model data

    using ModAstroConstants::M_Sun;
    using ModAstroConstants::L_Sun;
    using ModAstroConstants::R_Sun;
    using ModAstroConstants::four_pi;
    using ModAstroConstants::sigma;

    char y_n;

    //Write introductory information for user
    cout << "               StatStar is designed to build a ZAMS star\n\n\n"
        <<  "               Details of the code are described in:\n"
        <<  "                  An Introduction to Modern Astrophysics\n"
        <<  "                  Bradley W. Carroll and Dale A. Ostlie\n"
        <<  "                     Second Edition, Addison Wesley\n"
        <<  "                     copyright 2007\n" << endl;

    cout << "               The user will be asked to enter the following quantities:\n"
        <<  "                  Mass of the star       (in solar units)\n"
        <<  "                  Luminosity of the star (in solar units)\n"
        <<  "                  Effective Temperature  (in K)\n"
        <<  "                  Hydrogen mass fraction (X)\n"
        <<  "                  Metals mass fraction   (Z)\n" << endl;

    bool All_parameters = true;
    bool Fix_parameters = true;
    bool io_ok;

    while (All_parameters) {
        if (all_new == 'Y' || all_new == 'y' || all_new == 'A' || all_new == 'a') {
            do {
                cout << "Enter the mass (in solar units)       :  Msolar = ";
                cin  >> Msolar;
                if (cin.fail() || Msolar <= 0) {
                    io_ok = false;
                    cout << "Invalid value entered - please try again" << endl;
                    if (cin.fail()){
                        cin.clear();
                        char c;
                        cin >> c;}
                }
                else io_ok = true;
            } while (!io_ok);

            do {
                cout << "Enter the effective temperature (in K):  Teff   = ";
                cin  >> Teff;
                if (cin.fail() || Teff <= 0) {
                    io_ok = false;
                    cout << "Invalid value entered - please try again" << endl;
                    if (cin.fail()){
                        cin.clear();
                        char c;
                        cin >> c;}
                }
                else io_ok = true;
            } while (!io_ok);

            do {
                cout << "Enter the luminosity (in solar units) :  Lsolar = ";
                cin  >> Lsolar;
                if (cin.fail() || Lsolar <= 0) {
                    io_ok = false;
                    cout << "Invalid value entered - please try again" << endl;
                    if (cin.fail()){
                        cin.clear();
                        char c;
                        cin >> c;}
                }
                else io_ok = true;
            } while (!io_ok);

            do {
                do {
                    cout << "Enter the mass fraction of hydrogen   :  X      = ";
                    cin  >> X;
                    if (cin.fail() || X < 0 || X > 1) {
                        io_ok = false;
                        if (X < 0 || X > 1) cout << "\n0 <= X <= 1 is required";
                        cout << "\nInvalid value entered - please try again\n" << endl;
                        if (cin.fail()){
                            cin.clear();
                            char c;
                            cin >> c;}
                    }
                    else io_ok = true;
                } while (!io_ok);

                do {
                    cout << "Enter the mass fraction of metals     :  Z      = ";
                    cin  >> Z;
                    if (cin.fail() || Z < 0 || Z > 1) {
                        io_ok = false;
                        if (Z < 0 || Z > 1) cout << "\n0 <= Z <= 1 is required";
                        cout << "\nInvalid value entered - please try again\n" << endl;
                        if (cin.fail()){
                            cin.clear();
                            char c;
                            cin >> c;}
                    }
                    else io_ok = true;
                } while (!io_ok);

                Y = Helium(X, Z);
                if (Y < 0 || Y > 1) {
                    cout << "Note that 0 <= X + Z <= 1 is required\n"
                        <<  "  Please reenter composition values" << endl;
                io_ok = false;}
            } while (!io_ok);
        }

        //Compute SI values
        Ms = Msolar*M_Sun;
        Ls = Lsolar*L_Sun;
        Rs = sqrt(Ls/(four_pi*sigma*pow(Teff, 4)));     //Eq. (3.17)
        Rsolar = Rs/R_Sun;                              //Solar radius as SI value

        //Allow the user the opportunity to change values as needed
        while (Fix_parameters) {
            cout << fixed << showpoint << right << setprecision(5) << "\n"
                << "Your model parameters are:               "
                << "M      = " << setw(11) << Msolar << " solar" << "\n"
                << "                                         "
                << "Teff   = " << setw(11) << Teff   << " K    " << "\n"
                << "                                         "
                << "L      = " << setw(11) << Lsolar << " solar" << "\n"
                << "                                         "
                << "R      = " << setw(11) << Rsolar << " solar" << "\n"
                << "                                         "
                << "X      = " << setw(11) << X << "\n"
                << "                                         "
                << "Y      = " << setw(11) << Y << "\n"
                << "                                         "
                << "Z      = " << setw(11) << Z << endl;

            cout << "\n" << "Are these values ok (y/n)? ";
            do {
                cin >> y_n;
                if (cin.fail() || !(y_n == 'Y' || y_n == 'y' || y_n == 'N' || y_n == 'n')){
                    io_ok = false;
                    cout << "Please answer Yes (y) or No (n):  ";
                    if (cin.fail()){
                        cin.clear();
                        char c;
                        cin >> c;}
                }
                else io_ok = true;
            } while (!io_ok);
            
            if (y_n == 'Y' || y_n == 'y') {
                All_parameters = false;
                Fix_parameters = false;}
            else {
                all_new = 'N';
                Change_Model(Msolar, Lsolar, Rsolar, Ms, Ls, Rs, Teff, X, Y, Z, all_new);
                if (all_new == 'E' || all_new == 'e') {
                    All_parameters = false;
                    Fix_parameters = false;}
                if (all_new == 'A' || all_new == 'a') Fix_parameters = false;}
        }
    }
    return;
}



void Change_Model(double& Msolar, double& Lsolar, double& Rsolar, double& Ms, double& Ls, 
                   double& Rs, double& Teff, double& X, double& Y, double& Z, char& all_new)
{
//      General Description:
//      ====================
//          Get updated model input data

    using ModAstroConstants::M_Sun;
    using ModAstroConstants::L_Sun;
    using ModAstroConstants::R_Sun;
    using ModAstroConstants::four_pi;
    using ModAstroConstants::sigma;

    char y_n;
    bool io_ok;

    if (all_new == 'Y' || all_new == 'y') {
        cout << "\n" << "Would you like to run another model?\n"
            << "Your previous results will be overwritten in the output file. (y/n): ";
        do {
            cin >> y_n;
            if (cin.fail() || !(y_n == 'Y' || y_n == 'y' || y_n == 'N' || y_n == 'n')){
                io_ok = false;
                cout << "Please answer Yes (y) or No (n):  ";
                if (cin.fail()){
                    cin.clear();
                    char c;
                    cin >> c;}
            }
            else io_ok = true;
        } while (!io_ok);

        if (y_n == 'Y' || y_n == 'y') {
            cout << fixed << showpoint << right << setprecision(5)
                << "\nWhich variable would you like to change?\n"
                << "     M = Mass                          Current value = " << setw(11) << Msolar << " solar\n"
                << "     T = effective Temperature         Current value = " << setw(11) << Teff   << " K    \n"
                << "     L = Luminosity                    Current value = " << setw(11) << Lsolar << " solar\n"
                << "     X = hydrogen mass mraction (X)    Current value = " << setw(11) << X << "\n"
                << "     Z = metals mass fraction (Z)      Current value = " << setw(11) << Z << "\n"
                << "     A = select an All new set of model parameters\n"
                << "     E = Exit the calculation\n" << endl;
            cout << "Select a letter: ";

            do {
                cin >> all_new;
                if (cin.fail() || 
                   !(all_new == 'M' || all_new == 'm' ||
                     all_new == 'T' || all_new == 't' ||
                     all_new == 'L' || all_new == 'l' ||
                     all_new == 'X' || all_new == 'x' ||
                     all_new == 'Z' || all_new == 'z' ||
                     all_new == 'A' || all_new == 'a' ||
                     all_new == 'E' || all_new == 'e')) {
                         
                    io_ok = false;
                    cout << "Please respond with one of the options listed:  ";
                    if (cin.fail()){
                        cin.clear();
                        char c;
                        cin >> c;}
                }
                else io_ok = true;
            } while (!io_ok);
        }
    }
    else {
        y_n = 'Y';
        cout << "\nWhich variable would you like to change?\n"
            << "     M = Mass\n"
            << "     T = effective Temperature\n"
            << "     L = Luminosity\n"
            << "     X = hydrogen mass fraction (X)\n"
            << "     Z = metals mass fraction (Z)\n"
            << "     A = select an All new set of model parameters\n"
            << "     E = Exit the calculation\n" << endl;
        cout << "Select a letter: ";

        do {
            cin >> all_new;
            if (cin.fail() || 
                !(all_new == 'M' || all_new == 'm' ||
                  all_new == 'T' || all_new == 't' ||
                  all_new == 'L' || all_new == 'l' ||
                  all_new == 'X' || all_new == 'x' ||
                  all_new == 'Z' || all_new == 'z' ||
                  all_new == 'A' || all_new == 'a' ||
                  all_new == 'E' || all_new == 'e')) {

                io_ok = false;
                cout << "Please respond with one of the options listed:  ";
                if (cin.fail()){
                    cin.clear();
                    char c;
                    cin >> c;}
            }
            else io_ok = true;
            } while (!io_ok);

        if (all_new == 'E' || all_new == 'e') y_n = 'N';
    }

    if (y_n == 'Y' || y_n == 'y') {
        switch (all_new) {
            case 'M':
            case 'm':
                do {
                    cout << "Enter the mass (in solar units)       :  Msolar = ";
                    cin >> Msolar;
                    if (cin.fail() || Msolar <= 0) {
                        io_ok = false;
                        cout << "Invalid value entered - please try again" << endl;
                        if (cin.fail()){
                            cin.clear();
                            char c;
                            cin >> c;}
                    }
                    else {
                        io_ok = true;
                        Ms = Msolar*M_Sun;}
                } while (!io_ok);
                all_new = 'n';
                break;
            
            case 'T':
            case 't':
                do {
                    cout << "Enter the effective temperature (in K):  Teff   = ";
                    cin >> Teff;
                    if (cin.fail() || Teff <= 0) {
                        io_ok = false;
                        cout << "Invalid value entered - please try again" << endl;
                        if (cin.fail()){
                            cin.clear();
                            char c;
                            cin >> c;}
                    }
                    else {
                        io_ok = true;
                        Rs = sqrt(Ls/(four_pi*sigma*pow(Teff, 4)));     //Eq. (3.17)
                        Rsolar = Rs/R_Sun;}                             //Solar radius from SI value
                } while (!io_ok);
                all_new = 'n';
                break;

            case 'L':
            case 'l':
                do {
                    cout << "Enter the luminosity (in solar units) :  Lsolar = ";
                    cin >> Lsolar;
                    if (cin.fail() || Lsolar <= 0) {
                        io_ok = false;
                        cout << "Invalid value entered - please try again" << endl;
                        if (cin.fail()){
                            cin.clear();
                            char c;
                            cin >> c;}
                    }
                    else {
                        io_ok = true;
                        Ls = Lsolar*L_Sun;
                        Rs = sqrt(Ls/(four_pi*sigma*pow(Teff, 4)));     //Eq. (3.17)
                        Rsolar = Rs/R_Sun;}                             //Solar radius from SI value
                } while (!io_ok);
                all_new = 'n';
                break;

            case 'X':
            case 'x':
                do {
                    cout << "Enter the mass fraction of hydrogen   :  X      = ";
                    cin >> X;
                    if (cin.fail() || X < 0 || X > 1) {
                        io_ok = false;
                        if (X < 0 || X > 1) cout << "\n0 <= X <= 1 is required" << endl;
                        cout << "Invalid value entered - please try again\n" << endl;
                        if (cin.fail()){
                            cin.clear();
                            char c;
                            cin >> c;}
                    } else io_ok = true;
                } while (!io_ok);

                Y = Helium(X, Z);
                if (Y < 0 || Y > 1) {
                    do {
                        cout << "Note that 0 <= X + Z <= 1 is required" << endl;
                        do {
                            cout << "  Please try again:                      X      = ";
                            cin >> X;
                            if (cin.fail() || X < 0 || X > 1) {
                                io_ok = false;
                                if (X < 0 || X > 1) cout << "\n0 <= X <= 1 is required" << endl;
                                cout << "Invalid value entered\n" << endl;
                                if (cin.fail()){
                                    cin.clear();
                                    char c;
                                    cin >> c;}
                            } else io_ok = true;
                        } while (!io_ok);
                        Y = Helium(X, Z);
                    } while (Y < 0 || Y > 1);
                }
                all_new = 'n';
                break;

            case 'Z':
            case 'z':
                do {
                    cout << "Enter the mass fraction of metals     :  Z      = ";
                    cin >> Z;
                    if (cin.fail() || Z < 0 || Z > 1) {
                        io_ok = false;
                        if (Z < 0 || Z > 1) cout << "\n0 <= Z <= 1 is required" << endl;
                        cout << "Invalid value entered - please try again\n" << endl;
                        if (cin.fail()){
                            cin.clear();
                            char c;
                            cin >> c;}
                    } else io_ok = true;
                } while (!io_ok);

                Y = Helium(X, Z);
                if (Y < 0 || Y > 1) {
                    do {
                        cout << "Note that 0 <= X + Z <= 1 is required" << endl;
                        do {
                            cout << "  Please try again:                      Z      = ";
                            cin >> X;
                            if (cin.fail() || Z < 0 || Z > 1) {
                                io_ok = false;
                                if (Z < 0 || Z > 1) cout << "\n0 <= Z <= 1 is required" << endl;
                                cout << "Invalid value entered\n" << endl;
                                if (cin.fail()){
                                    cin.clear();
                                    char c;
                                    cin >> c;}
                            } else io_ok = true;
                        } while (!io_ok);
                        Y = Helium(X, Z);
                    } while (Y < 0 || Y > 1);
                }
                all_new = 'n';
                break;

            case 'E':
            case 'e':
                y_n = 'N';
                break;

            case 'A':
            case 'a':
                y_n = 'Y';
                break;

            default:
                all_new = 'y';
        }
    }
    else all_new = 'E';             //Exit calculations

    return;
}



void Too_Many_Zones(int i, double Msolar, double Ms, double M_r, double Lsolar, double Ls, double L_r, 
                    double r, double Rs, double Rsolar, double Teff, double X, double Y, double Z, 
                    double P_0, double T_0, double rho_0, double kappa_0, double epsilon_0, char rc_flag)
{
//      General Description:
//      ====================
//          Tell user that the maximum number of zones has been exceeded

    cout << fixed << showpoint << right << setprecision(6) << "\n"
        << "The maximum number of zones has been exceeded for this model - Sorry\n"
        << "The conditions at the time the model was terminated were: \n"
        << "    Surface Conditions:            Last Zone Calculated:\n"
        << "    -------------------            ---------------------\n"
        << "    M    = " << setw(13) << Msolar << " solar     M_r/Ms  = " << setw(13) << M_r/Ms << "\n"
        << "    Teff = " << setw(13) << Teff   << " K         L_r/LS  = " << setw(13) << L_r/Ls << "\n"
        << "    L    = " << setw(13) << Lsolar << " solar     r/Rs    = " << setw(13) << r/Rs << "\n"
        << "    R    = " << setw(13) << Rsolar << " solar     P       = " 
        << scientific << setprecision(5) 
        <<                                                                   setw(13) << P_0 << " N/m^2\n"
        << fixed << setprecision(6)
        << "    X    = " << setw(13) << X      << "           T       = "
        << scientific << setprecision(5)
        <<                                                                   setw(13) << T_0 << " K\n"
        << fixed << setprecision(6)
        << "    Y    = " << setw(13) << Y      << "           rho     = "
        << scientific << setprecision(5)
        <<                                                                   setw(13) << rho_0 << " kg/m^3\n"
        << fixed << setprecision(6)
        << "    Z    = " << setw(13) << Z      << "           kappa   = "
        << scientific << setprecision(5)
        <<                                                                   setw(13) << kappa_0 << " m^2/kg\n"
        <<                "                                   epsilon = " << setw(13) << epsilon_0 << " W/kg"
        << endl;

    if (rc_flag == 'r') cout << "                                  The core is RADIATIVE" << endl;
    else                cout << "                                  The core is CONVECTIVE" << endl;

    return;
}

void Final_Results(int i, double Msolar, double Ms, double M_r, double Lsolar, double Ls, double L_r, 
                   double r, double Rs, double Rsolar, double Teff, double X, double Y, double Z, 
                   double P, double T, double rho, double kappa, double epsilon, 
                   double P_0, double T_0, double rho_0, double kappa_0, double epsilon_0, char rc_flag)
{
//      General Description:
//      ====================
//          Tell the user the conditions at the surface and core of the completed model


    cout << fixed << showpoint << right << setprecision(6) << "\n\n"
        << "*********************THE INTEGRATION HAS BEEN COMPLETED*********************\n"
        << "    Surface Conditions:            Last Zone Calculated:\n"
        << "    -------------------            ---------------------\n"
        << "    M    = " << setw(13) << Msolar << " solar     M_r/Ms  = " << setw(13) << M_r/Ms << "\n"
        << "    Teff = " << setw(13) << Teff   << " K         L_r/LS  = " << setw(13) << L_r/Ls << "\n"
        << "    L    = " << setw(13) << Lsolar << " solar     r/Rs    = " << setw(13) << r/Rs << "\n"
        << "    R    = " << setw(13) << Rsolar << " solar     P       = " 
        << scientific << setprecision(5) 
        <<                                                                   setw(13) << P_0 << " N/m^2\n"
        << fixed << setprecision(6)
        << "    X    = " << setw(13) << X      << "           T       = "
        << scientific << setprecision(5)
        <<                                                                   setw(13) << T_0 << " K\n"
        << fixed << setprecision(6)
        << "    Y    = " << setw(13) << Y      << "           rho     = "
        << scientific << setprecision(5)
        <<                                                                   setw(13) << rho_0 << " kg/m^3\n"
        << fixed << setprecision(6)
        << "    Z    = " << setw(13) << Z      << "           kappa   = "
        << scientific << setprecision(5)
        <<                                                                   setw(13) << kappa_0 << " m^2/kg\n"
        <<                "                                   epsilon = " << setw(13) << epsilon_0 << " W/kg"
        << endl;

    if (rc_flag == 'r') cout << "                                   The core is RADIATIVE" << endl;
    else                cout << "                                   The core is CONVECTIVE" << endl;

    cout << scientific << showpoint << right << setprecision(5)
        << "\nFor your information, the conditions in the last zone above the core are:\n"
        << "                                   P       = " << setw(13) << P       << " N/m^2\n"
        << "                                   T       = " << setw(13) << T       << " K\n"
        << "                                   rho     = " << setw(13) << rho     << " kg/m^3\n"
        << "                                   kappa   = " << setw(13) << kappa   << " m^2/kg\n"
        << "                                   epsilon = " << setw(13) << epsilon << " W/kg" << endl;

    cout << fixed << noshowpoint 
        << "\nThe number of mass shells in this model: " << setw(5) << i
        << "\nThe details of the model are available in ZAMSmodel.txt" << endl;
}
