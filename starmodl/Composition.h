/* Composition

   General Description:
   ====================
       This header file contains inlined function definitions for information about the composition of the gas

---------------------------------------------------------------------*/

#ifndef COMPOSITION_H
#define COMPOSITION_H

using namespace std;

//      General Description:
//      ====================
//          Calculate the mean molecular weight of the gas

inline double Mean_Molecular_Weight(double X, double Y, double Z)
{
    return (1/(2*X + 3*Y/4 + Z/2));         //Assume complete ionization, Eq. (10.16)
}

//--------------------------------------------------------------------

//      General Description:
//      ====================
//          Calculate the amount of Helium-4 in the mixture

inline double Helium(double X, double Z)
{
    return (1 - X - Z);
}

//--------------------------------------------------------------------

//      General Description:
//      ====================
//          Calculate the mass fraction of C, N, and O in the mixture

inline double CNO(double Z)
{
    return (Z/2);
}

#endif
