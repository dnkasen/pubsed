import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys


def run_test(pdf="",runcommand=""):

    #-------------------------------------------
    """ Function to run a test of the sedona code

        Args:
            pdf: pointer to a pdf file to output data to
            runcommand: string giving the command to run,
                        e.g., "mpirun -np 6 ./sedona6"

        Returns:
            an integer ("failure") specifying success or failure of test
            failure == 0 if success
            failure != 0 if failed (with the number being some code for what failed)

    """
    #-------------------------------------------

    testname = "1D toy type1a supernova spectrum - exp"

    #-------------------------------------------
    # clean up any old results and run the code
    #-------------------------------------------
    if (runcommand != ""):
    	os.system("rm spectrum_* plt_* integrated_quantities.dat")
    	os.system(runcommand)

    #-------------------------------------------
    # compare the output
    #-------------------------------------------
    failure = 0
    plt.clf()

    # plot the spectrum
    nu,Lnu = np.loadtxt('spectrum_5.dat',unpack=True,skiprows=1,usecols=[0,1])
    lam = 3e10/nu*1e8
    Llam = Lnu*nu/lam
    plt.plot(lam,Llam,color='k',lw=2)

    # plot the comparison spectrum
    nu,Lnu = np.loadtxt('reference_spectrum.dat',unpack=True,skiprows=1,usecols=[0,1])
    lam = 3e10/nu*1e8
    Llam_ref = Lnu*nu/lam
    plt.plot(lam,Llam_ref,color='r',lw=2)

    plt.xlim(1000,10000)
    plt.xlabel('wavelength (angstroms)')
    plt.ylabel('specific luminoisty (ergs/s/angstrom)')
    plt.legend(['sedona output','reference'])
    plt.title(testname)

    # write code here to do a numerical comparison
    # of output data to standard reference and
    # determine success or failure
    # set failure to some number != 0 if something
    # is wrong

    use = (Llam_ref > max(Llam_ref)*0.1)
    max_err,mean_err = get_error(Llam,Llam_ref,use=use)
    if (mean_err > 0.1): failure = 1

    # add plot to pdf file (or show on screen)
    if (pdf != ''): pdf.savefig()
    else:
        plt.ion()
        plt.show()
        j = get_input('Press any key to continue >')

    # this should return !=0 if failed
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
