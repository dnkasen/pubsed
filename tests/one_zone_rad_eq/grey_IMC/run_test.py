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
    	os.system("rm spectrum_* plt_* integrated_quantities.dat")
    	os.system(runcommand)

    failure = 0
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
    gas_temp_initial = 5.5e8

    arad = 4. * sb/ c
    pltfile_list = glob.glob('./plt_*dat')
    times = np.array([])
    egas = np.array([])

    npts = 25
    times = np.zeros(npts)
    egas  = np.zeros(npts)

    for i in range(0,npts):
        pltfile = 'plt_000' + str(i+1) + '.dat'
        if (i < 9):  pltfile = 'plt_0000' + str(i+1) + '.dat'
        times[i] = float((open(pltfile, 'r').readline()).split()[3])
        data = np.loadtxt(pltfile,skiprows=2)
        rho = data[1]
        gas_temp = data[3]
        egas[i] = 1./(gamma - 1.) * rho / (mean_particle_mass * mp) * k * gas_temp


    #------------------------------------------
    #compare output 
    #------------------------------------------

    plt.clf()

    plt.plot(times,egas,'o',color='black')

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
    e_analytic = gas_internal_energy_analytic[:,0]

    plt.plot(physical_time_grid,e_analytic,'-',c = 'r')

    use = (times > 1e-9)
    max_err,mean_err = get_error(egas,e_analytic,x=times,x_comp=physical_time_grid,use = use)
    if (mean_err > 0.01): failure = 1
    if (max_err > 0.02): failure = 1


    plt.legend(['sedona','analytic'])
    plt.title('one zone radiative equilibrium test')
    plt.xlabel('time (s)')
    plt.ylabel('energy density (erg / cm^3)')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(1.e-12,1.e-6)
    #plt.ylim(1.e-8,1.e-1)

    plt.show()
    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = raw_input("Press Return to continue>")


    return failure

#----------------------------------------------
# helper function for calculating the error
# between two numpy arrays
#----------------------------------------------

def get_error(a,b,x=[],x_comp=[],use=[]):

    #-------------------------------------------

    """ Function to calculate the error between two arrays

        Args:
        a: numpy array of result
        b: numpy array of comparison
        use: an array of 0's and 1's telling which element
             in the arrays to include
        x: optional array of x values to go along with a
        x_comp: optional array of x values to go along with b
        (if x and x_comp are set, will interpolate b values to x spacing)

        Returns:
            returns max_error, mean_error in percentages

        Example:
            say you have an array y that is a function of x
            you wnat to see how much it deviates from a reference array y_comp
            but only for values where x > 0.5. Use

            max_error, mean_error = get_error(y,y_comp,use=(x > 0.5))

    """
    #-------------------------------------------------

    # result array
    y = a
    # compare array
    y_comp = b

    # interpolate comparison if wanted
    if (len(x) != 0 and len(x_comp !=0)):
        y_comp = np.interp(x,x_comp,y_comp)

    # cut the array length if wanted
    if (len(use) > 0):
        y = y[use]
        y_comp = y_comp[use]
    err = abs(y - y_comp)

    max_err = max(err/y_comp)
    mean_err = np.mean(err)/np.mean(y_comp)

    return max_err,mean_err


#-----------------------------------------
# little function to just plot up and
# compare results. Assumes code has
# already been run and output files
# are present
#----------------------------------------
if __name__=='__main__':

    # Default to Python 3's input()
    get_input = input
    # If this is Python 2, use raw_input()
    if sys.version_info[:2] <= (2, 7):
        get_input = raw_input

    status = run_test('')
    if (status == 0):
        print ('SUCCESS')
    else:
        print ('FAILURE, code = ' + str(status))