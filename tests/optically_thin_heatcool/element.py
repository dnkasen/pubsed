#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python
import constants as const
import numpy as np
from math import pow, exp, sqrt
import atomic_level as level
import atomic_line as line
import photoion_cs

from scipy.integrate import simps, trapz

#################

class Element:

    def __init__(self, Z):

        self.atomic_number = Z
        self.num_levels = 0
        self.num_lines = 0
        self.levels = []
        self.lines = []
        self.total_ionization_energy_ev = 13.60 * pow(Z,2) # energy to remove the last electron
        self.max_ionization_cs = 6.30e-18 / pow(Z,2) # for highest ionization state


    def set_total_ionization_energy_ev(self, energy):
        self.total_ionization_energy_ev = energy
        
    def add_level(self,i,n,E,g): 
 
        new_level = level.AtomicLevel(i,n,E,g)
        self.levels.append(new_level)
        self.num_levels += 1
        self.levels[-1].ionization_frequency_threshold = (self.total_ionization_energy_ev - self.levels[-1].excitation_energy_ev) * const.ev_to_ergs/const.h 

    def add_line(self,i,li,ui,A_ul):
        self.lines.append(line.AtomicLine(i,li,ui,A_ul))
        self.num_lines += 1
        self.lines[-1].nu = (self.levels[ui].excitation_energy_ev - self.levels[li].excitation_energy_ev) * const.ev_to_ergs/const.h

    def approx_photoion_cs(self,nu,i,n_index):
        nu_nd = nu / self.levels[i].ionization_frequency_threshold 

        if (nu_nd == 1.):
            return n_index * self.max_ionization_cs
        else:
            if (nu_nd > 1.):
                return n_index * self.max_ionization_cs * pow(nu_nd,-3.)
            else:
                return 0

    def exact_photoion_cs_groundstate(self,nu):
        return photoion_cs.exact_a_nu_1(nu)


#Again, just using hydrogenic approximation to set effective "n" for each level. And this will need to be expanded when you keep track of ionization states of elements besides hydrogen
    def compute_alpha(self, i, gas_temp,bbobj):

        #Maxwell-Boltzmann constants

        n = self.levels[i].n_quantum

        v_MB = sqrt(2.*const.k*gas_temp/const.m_e)
        MB_A = 4./sqrt(np.pi)*pow(v_MB,-3.)
        milne_prefac = pow(const.h/const.c/const.m_e,2.) 

        nu_t  = (self.total_ionization_energy_ev - self.levels[i].excitation_energy_ev) * const.ev_to_ergs/const.h

        integrand = np.zeros(len(bbobj.frequencies))

        for index in range(len(bbobj.frequencies)):

            nu  = bbobj.frequencies[index]

            if (nu >= nu_t):

                if (index == 0):
                    S   = self.exact_photoion_cs_groundstate(bbobj.frequencies[index])
                else:
                    S   = self.approx_photoion_cs(bbobj.frequencies[index],i,n)
                
                u_e  = sqrt(2.*const.h*(nu - nu_t)/const.m_e)

                f_u = MB_A * u_e * u_e * exp(-1. * const.m_e * u_e * u_e/(2. * const.k * gas_temp)) 

                sigma = milne_prefac * nu * nu / (u_e * u_e) * S

                integrand[index]  = u_e*sigma*f_u * sqrt(2. *const.h/ const.m_e) * 1./2. * pow(nu - nu_t,-1./2.)

            else:
                
                integrand[index]  = 0.

        # integrate

        #integral = trapz(integrand, bbobj.frequencies)
        integral = simps(integrand, bbobj.frequencies)

        return self.levels[i].statistical_weight* integral       

                                       
        







