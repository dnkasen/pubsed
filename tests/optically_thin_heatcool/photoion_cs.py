#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import constants as const
import numpy as np
from math import pow, sqrt, exp, atan


#################

# Cross-section for hydrogen photoionization from ground state, including quantum corrections

nu_threshold_hydrogen_1 = 3.29e15
sigma_threshold_hydrogen = 6.30e-18

nu_threshold_hydrogen_2 = (13.59808 - 10.198838) * const.ev_to_ergs / const.h
nu_threshold_hydrogen_3 = (13.59808 - 12.087516) * const.ev_to_ergs / const.h

def epsilon(nu):
    return sqrt(nu - 1)

def exact_a_nu_1(nu):
    nu_nd = nu / nu_threshold_hydrogen_1
    if (nu_nd == 1.):
        return sigma_threshold_hydrogen
    else:
        if (nu_nd > 1.):
            return sigma_threshold_hydrogen * pow(nu_nd,-4) * exp(4 - ((4 * atan(epsilon(nu_nd)))/epsilon(nu_nd)))/(1 - exp(-2 * np.pi/epsilon(nu_nd)))
        else:
            return 0.

def approx_a_nu(nu,n_index):
    if (n_index == 1):
        nu_nd = nu / nu_threshold_hydrogen_1 
    if (n_index == 2):
        nu_nd = nu / nu_threshold_hydrogen_2
    if (n_index == 3):
        nu_nd = nu / nu_threshold_hydrogen_3

    if (nu_nd >= 1.):
        return n_index * sigma_threshold_hydrogen * pow(nu_nd,-3.)
    else:
        return 0

def approx_a_nu_1(nu):
    nu_nd = nu /nu_threshold_hydrogen_1
    if (nu_nd == 1.):
        return 1. * sigma_threshold_hydrogen
    else:
        if (nu_nd > 1.):
            return 1. * sigma_threshold_hydrogen * pow(nu_nd,-3.)
        else:
            return 0


def approx_a_nu_2(nu):
    nu_nd = nu / nu_threshold_hydrogen_2 
    if (nu_nd == 1.):
        return 2. * sigma_threshold_hydrogen
    else:
        if (nu_nd > 1.):
            return 2. * sigma_threshold_hydrogen * pow(nu_nd,-3.)
        else:
            return 0

def approx_a_nu_3(nu):
    nu_nd = nu / nu_threshold_hydrogen_3
    if (nu_nd == 1.):
        return 3. * sigma_threshold_hydrogen
    else:
        if (nu_nd > 1.):
            return 3. * sigma_threshold_hydrogen * pow(nu_nd,-3.)
        else:
            return 0

#################

# Cross-section for photoionization of first electron of helium from ground state, based on Verner et al (1996) fit 

sigma_0 = 9.492e-16 # cm^2
y_0 = 0.4434
E_0 = 2.181e-11 # ergs
y_1 = 2.136
Pexp = 3.188
y_w = 2.039 
y_a = 1.469 

def x_E(energy):
    return energy/E_0 - y_0
    
def y_E(xa):
    return sqrt(pow(xa,2.) + pow(y_1,2.))
    
def sigma_F(nu):    
    return (pow(x_E(const.h * nu) - 1.,2.) + pow(y_w,2.)) * pow(y_E(x_E(const.h * nu)), 0.5 * Pexp - 5.5) * pow(1. + sqrt(y_E(x_E(const.h * nu))/y_a), -Pexp)
    
def anu_helium(nu):
    if (nu < nu_threshold_helium):
        return 0.
    else:
        return sigma_0 * sigma_F(nu) 

###################

# Cross-section of second ionization of helium from ground state. This is just like hydrogen, but appropriately rescaled for new value of Z

def anu_heplus(nu):
    nu_nd = nu / (nu_threshold_hydrogen * 4.)
    if (nu_nd == 1.):
        return sigma_threshold_hydrogen / 4.
    else:
        if (nu_nd > 1.):
            return sigma_threshold_hydrogen/4. * pow(nu_nd,-4) * exp(4 - ((4 * atan(epsilon(nu_nd)))/epsilon(nu_nd)))/(1 - exp(-2 * np.pi/epsilon(nu_nd)))
        else:
            return 0.
