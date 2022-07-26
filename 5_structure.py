# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# Time-stamp: <2016-05-11 10:21:58 (jmiller)>

# This is a simple script to calculate and plot the pressure in a star
# as a function of radius.

# To run use
#python3 structure.py

#This should produce a pdf plot. You can play with the constants at the
#top of the program to change the polytropic index, the initial
#pressure, or the radial distance for the forward Euler step.

# Imports
# ------------------------------------------------
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
# ------------------------------------------------

# Constants
# ------------------------------------------------
K=1.0 # polytrope K
n=3.0 # polytropic index
exponent = ((n+1)/n)
G=1.0 # Newton's constant
# ------------------------------------------------


def get_rho(P):
    """The equation of state gives ensity a function of pressure.
    We assume a polytropic equation of state:
    P = K*rho**((n+1)/n)
    which we invert"""
    out = P/K
    rho = np.abs(out)**(1/exponent)
    return rho

def rhs(r,v):
    """Right-hand side
    takes in vector (M,P)
    and returns dM/dr, dP/dr
    also requires radius r"""
    M,P = v
    rho = get_rho(P)
    dmdr=4*np.pi*r*r*rho
    dpdr = - G*M*rho/(r*r)
    return np.array([dmdr,dpdr])

def forward_euler_step(r,v,dr):
    """Uses the stellar structure equations
    to integrate outward from the center of the star"""
    dv = dr*rhs(r,v)
    return v + dv

def integrate_test(p0,dr):
    """Gives us the final radius to integrate to."""
    M0 = 0
    r0 = dr
    P0 = p0
    v = np.array([M0,P0])
    r = r0
    nr = 0
    while v[1] > 0:
        v = forward_euler_step(r,v,dr)
        r +=  dr
        nr += 1
    return r,nr

def integrate_final(p0,dr,nr):
    """Once we know the final radius, we integrate to it
    and store our data."""
    M0 = 0
    r0 = dr
    P0 = p0
    v0 = np.array([M0,P0])
    r = np.arange(dr,(nr+1)*dr,dr)
    v = np.empty((nr,2),dtype=float)
    v[0,...] = v0
    for i in range(1,nr):
        v[i,...] = forward_euler_step(r[i-1],v[i-1],dr)
    return r,v

def plot_profile(r,v):
    m = v[...,0]
    p = v[...,1]
    rho = get_rho(p)
    plt.plot(r,m,lw=3)
    plt.plot(r,p,lw=3)
    plt.plot(r,rho,lw=3)
    plt.xlabel('radius',fontsize=20)
    plt.legend(['mass','pressure','density'],
               fontsize=20,
               loc = 'lower right')
    plt.savefig('star_profile.pdf',
                bbox_inches='tight')
    
if __name__ == "__main__":
    p0 = 1.4
    dr = 0.001
    rmax,nr = integrate_test(p0,dr)
    r,v = integrate_final(p0,dr,nr)
    plot_profile(r,v)
