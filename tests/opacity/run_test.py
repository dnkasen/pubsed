import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys




def run_test(pdf="",runcommand=""):

    ###########################################
    # clean up old results and run the code
    ###########################################
    if (runcommand != ""):
    	os.system("rm spectrum_* plt_* integrated_quantities.dat")
    	os.system(runcommand)

    ###########################################
    # compare the output
    ###########################################
    plt.clf()

    # plot output
    pltfile = h5py.File('plt_00001.h5', 'r')
    nu = np.array(pltfile["nu"])
    opac = np.array(pltfile["zonedata/0/opacity"])
    lam = 3e18/nu
    plt.plot(lam,opac,lw=3,color='k')
    pltfile.close()

    # plot reference
    pltfile = h5py.File('reference.h5', 'r')
    nu_r = np.array(pltfile["nu"])
    opac_r = np.array(pltfile["zonedata/0/opacity"])
    lam_r = 3e18/nu_r
    plt.plot(lam_r,opac_r,lw=3,color='r')

    plt.yscale('log')
    plt.xlabel('wavelength (angstroms)',size=12)
    plt.ylabel('Fe line expansion opacity (cm^2/g)',size=12)
    plt.title('opacity test; line expansion opacity')
    plt.legend(['Sedona','Reference'])

    use = ((nu > 3e14)*(nu < 2e15))
    max_err,mean_err = get_error(opac,opac_r,x=nu,x_comp=nu_r,use = use)
    #print(max_err,mean_err)

    if (pdf != ''): pdf.savefig()
    else:
        plt.show()

    # this should return !=0 if failed
    if (max_err > 0.05 or mean_err > 0.02):
        return 1
    return 0

#-------------------------------------------
# error calculator helper function
#-------------------------------------------

def get_error(a,b,x=[],x_comp=[],use=[]):

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

if __name__=='__main__':
    status = run_test('')
    if (status == 0):
        print ('SUCCESS')
    else:
        print ('FAILURE, code = ' + str(status))
