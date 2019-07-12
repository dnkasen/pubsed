#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import constants as const
import numpy as np
from math import pow, exp, log10, sqrt
from scipy.integrate import quad


class BlackbodyEmitter:


    def __init__(self):

        # params (default values)

        self.T_emit = 1.e4
        self.color_correction = 1.
        self.R_emit = 3.e13 
        self.longest_wavelength_angs = 1.e5
        self.shortest_wavelength_angs = 1.e1

        self.num_frequencies = 1600
        self.num_wavelengths = 1600

        self.reset_frequency_grid()
        self.reset_wavelength_grid()

#################

    def reset_frequency_grid(self):
        self.frequencies = np.logspace(np.log10(const.c / (self.longest_wavelength_angs * const.angs_to_cm) ),np.log10(const.c/(self.shortest_wavelength_angs * const.angs_to_cm)),self.num_frequencies)

    def reset_wavelength_grid(self):
        self.wavelengths = np.logspace(np.log10(self.shortest_wavelength_angs),np.log10(self.longest_wavelength_angs),self.num_wavelengths)

#        self.frequencies = np.linspace(const.c / (self.longest_wavelength_angs * const.angs_to_cm),const.c/(self.shortest_wavelength_angs * const.angs_to_cm),self.num_frequencies)

#################

    def PlanckFunctionFrequency(self,nu):
        try:
            return 2. * const.h * pow(nu,3.)/( pow(const.c,2) * (exp(const.h * nu /(const.k * self.T_emit)) - 1.) ) 
        except:
            if (const.h * nu /(const.k * self.T_emit) > 100.):
                return 0.
            else:
                print "error, check planck function"

    def PlanckFunctionWavelength(self,lam):
        try:
              return 2. * const.h * pow(const.c,2.)/( pow(lam,5) * (exp(const.h * const.c /(lam * const.k * self.T_emit)) - 1.) )
        except:
              if (const.h * const.c /(lam * const.k * self.T_emit) > 100.):
                  return 0.
              else:
                  print "Error, check Planck function"


    def dilution_factor(self,r):
        return 1./2. * (1. - sqrt(1. - pow(self.R_emit/r,2.) ))

    def EmittedJUnattenuated(self,nu,r):
        return self.PlanckFunctionFrequency(nu) * self.dilution_factor(r) * pow(self.color_correction,-4.)

    def PlanckFunctionFrequencyPhotons(self,nu):
        return self.PlanckFunctionFrequency(nu)/(const.h * nu)

    def SpecificPhotonEmissionRate(self,nu):
        return np.pi * self.PlanckFunctionFrequencyPhotons(nu) * 4. * np.pi * pow(self.R_emit,2) * pow(self.color_correction,-4.)

    def IntegratedPhotonEmissionRate(self,nu0):
        return quad(self.SpecificPhotonEmissionRate,nu0,self.frequencies[-1])[0]

    def TotalLuminosity(self):
        return 4. * np.pi * pow(self.R_emit,2.) * const.sb * pow(self.T_emit,4.) * pow(self.color_correction,-4.)

    #Mean energy for frequencies of nu0 and larger
    def MeanPhotonEnergy(self,nu0):
#        return (const.sb * pow(self.T_emit,4.) / np.pi) / quad(self.PlanckFunctionFrequencyPhotons,nu0,self.frequencies[-1])[0]
        return quad(self.PlanckFunctionFrequency,nu0,self.frequencies[-1])[0] / quad(self.PlanckFunctionFrequencyPhotons,nu0,self.frequencies[-1])[0]





