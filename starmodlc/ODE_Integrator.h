/* ODE_Integrator
 
    General Description:
    ====================
 
        ODE_Integrator contains the numerical integration routine used to integrate a set
        of n first order linear  differential equations that depend on one independent variable, x.
 
        This module accepts n ODEs with vectors of initial conditions.
        A user-defined function that returns the derivative of one of n specified functions is
        required to have the form:
 
            double Deriv(int i, double x, double y, bool& ok)

        where 
            (a)  all function derivative definitions are included in the single routine,
            (b)  i designates which function value is to be returned
            (c)  x is a scalar representing the independent variable
            (d)  y is a vector of dimension n representing the set of dependent variables at x
            (e)  ok is a logical flag.  If ok == .TRUE. the derivative was computed successfully.
 
---------------------------------------------------------------------*/

#ifndef ODE_INTEGRATOR_H
#define ODE_INTEGRATOR_H

void RK_4(int n, double x0, double h, double y0[], double y4[], double f0[], bool& ok,
          double X, double Z, double Pm, double Tm, int step_size_condition, char& rc_flag);

#endif
