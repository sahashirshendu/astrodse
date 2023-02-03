#include <iostream>
#include <iomanip>
#include "Stellar_Structure_Equations.h"

using namespace std;

void RK_4(int n, double x0, double h, double y0[], double y4[], double f0[], bool& ok,
          double X, double Z, double Pm, double Tm, int step_size_condition, char& rc_flag)
{
/*
        General Description:
        ====================
            This routine uses the fourth order Runge Kutta scheme to move the solution forward
            by one step size, h.
 
            For details of the method, see 
 
                Press, Teukolsky, Vetterling, and Flannery, "Numerical Recipes in C++:  The Art
                    of Scientific Computing", Second Edition, Cambridge University Press, Cambridge, 2002.
 
---------------------------------------------------------------------

        n               Number of ODEs
        x0, h           The independent variable and step size
        y0              The array of initial y values at the start of the interval
        y4              The array of results to 4th order at x + h
        f0              The array of first derivatives at the start of the interval
        ok              Reports if function evaluation was successful

        k1, k2, k3, k4  Temporary work arrays containing values at intermediate steps

!---------------------------------------------------------------------*/


    const int array_size = 4;
    double k1[array_size], k2[array_size], k3[array_size], k4[array_size], y0temp[array_size];

    if (n > array_size) {
        char xpause;
        cout << "Problem in RK_4 ... n exceeds size of temporary arrays" 
            "\nEnter any character to exit: "<< endl;
        cin >> xpause;
        exit(1);
    }

    int i;
    ok = true;

    //Calcualtion intermediate derivatives using the user-defined external function
    for (i = 0; i < n; i++){
        k1[i] = h*f0[i];}

    for (i = 0; i < n; i++){
        y0temp[i] = y0[i] + k1[i]/2;}

    for (i = 0; i < n; i++){
        k2[i] = h*Structure_Eqns(i, x0 + h/2, X, Z, Pm, Tm, step_size_condition, rc_flag, y0temp, ok);
        if (!ok) return;}

    for (i = 0; i < n; i++){
        y0temp[i] = y0[i] + k2[i]/2;}

    for (i = 0; i < n; i++){
        k3[i] = h*Structure_Eqns(i, x0 + h/2, X, Z, Pm, Tm, step_size_condition, rc_flag, y0temp, ok);
        if (!ok) return;}

    for (i = 0; i < n; i++){
        y0temp[i] = y0[i] + k3[i];}

    for (i = 0; i < n; i++){
        k4[i] = h*Structure_Eqns(i, x0 + h,   X, Z, Pm, Tm, step_size_condition, rc_flag, y0temp, ok);
        if (!ok) return;}

    //Compute the variables for the next shell using the 4th order Runge-Kutta formula
    for (i = 0; i < n; i++){
        y4[i] = y0[i] + k1[i]/6 + k2[i]/3 + k3[i]/3 + k4[i]/6;}

    return;
}