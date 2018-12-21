import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys
import glob
from scipy.integrate import odeint

def run_test(pdf="",runcommand=""):

    ###########################################
    # clean up old results and run the code
    ###########################################
    if (runcommand != ""):
    	os.system("rm *_spectrum_* plt_* integrated_quantities.dat")
    	os.system(runcommand)

    fail_flag = 0
    ###########################################
    # compare the output
    ###########################################
    plt.clf()


    h   = 6.6260755e-27    # planck's constant (ergs-s)
    c   = 2.99792458e10    # speed of light (cm/s)
    k   = 1.380658e-16     # boltzmann constant (ergs/K)
    sb  = 5.6704e-5        # stefan boltzman constant (ergs cm^-2 s^-1 K^-4)
    pi  = 3.14159          # just pi
    mp = 1.6726219e-24     # proton mass

    mean_particle_mass = 0.7
    gamma = 5./3.
    kappa = 0.4
    Tr = 3.4e6
    plt.ion()
    gas_temp_initial = 5.5e9

    arad = 4. * sb/ c
    pltfile_list = glob.glob('./plt_*dat')
    times = np.array([])
    egas = np.array([])
    for pltfile in pltfile_list:
        times = np.append(times,float((open(pltfile, 'r').readline()).split()[3]))
        data = np.loadtxt(pltfile,skiprows=2)
        rho = data[1]
        gas_temp = data[3]
        egas = np.append(egas,1./(gamma - 1.) * rho / (mean_particle_mass * mp) * k * gas_temp)

    #------------------------------------------
    #compare output 
    #------------------------------------------

    plt.clf()

    plt.plot(times,egas,'o',color='black')



#    plt.legend(['sedona','analytic blackbody'])
    plt.title('one zone radiative equilibrium test')
    plt.xlabel('time (s)')
    plt.ylabel('energy density (erg / cm^3)')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(1.e-12,1.e-6)
    #plt.ylim(1.e-8,1.e-1)


    # analytic solution
    # Using non-dimensionalized temperature theta and non-dimensional time tau
    #Define a function for the derivative
    def dthetadtau(theta,tau):
        return 1. - theta**4.
        #return (C_LIGHT * extinct * erad - 4 * S_BOLTZ * extinct * pow((gamfac - 1) * mu * M_PROTON * E/(rho * K_BOLTZ),4.))

    physical_time_grid = np.logspace(-16,-5,1024)
    characteristic_time = 1./(gamma - 1.) * rho / (mean_particle_mass * mp) * k /(arad * Tr**3. *  rho * kappa * c) 
    tau_grid = physical_time_grid / characteristic_time

    nondim_gastemp_analytic = odeint(dthetadtau,gas_temp_initial/Tr,tau_grid)
    gastemp_analytic = nondim_gastemp_analytic * Tr
    gas_internal_energy_analytic = 1./(gamma - 1.) * rho / (mean_particle_mass * mp) * k * gastemp_analytic

    plt.plot(physical_time_grid,gas_internal_energy_analytic,'-',c = 'r')

    plt.show()
    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = raw_input("Press Return to continue>")



    return fail_flag



if __name__=='__main__': run_test('')
