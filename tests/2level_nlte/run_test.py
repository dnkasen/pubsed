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

    testname = "my_test"

    #-------------------------------------------
    # clean up any old results and run the code
    #-------------------------------------------
    if (runcommand != ""):
    	os.system("rm *_spectrum_* plt_* integrated_quantities.dat")
    	os.system(runcommand)

    #-------------------------------------------
    # compare the output
    #-------------------------------------------
    failure = 0
    plt.clf()

    # write code here to read in output data and
    # plot it compared to standard reference

    fin = h5py.File('plt_00002.h5','r')
    rp = 3.0e13

    r    = np.array(fin['r'])
    rmid = r - 0.5*(r[1] - r[0])
    r   = rmid/rp
    rho = np.array(fin['rho'])
    nz = len(rho)

    n1 = np.zeros(nz)
    n2 = np.zeros(nz)
    for i in range(nz):
        name = 'zonedata/' + str(i)
        b = np.array(fin[name + '/Z_1/level_departure'])
        n = np.array(fin[name + '/Z_1/level_fraction'])
        n1[i] = n[0]
        n2[i] = n[1]
    plt.plot(r,n1/n2,'o',color='black')


    # analytic solution
    kb  = 1.380658e-16     # boltzmann constant (ergs/K)
    ev_to_ergs  = 1.60217646e-12
    T  = 5.0e4
    g1 = 2.0
    g2 = 8.0
    W = 0.5*(1 - (1 - (rp/rmid)**2)**0.5)
    eta = ev_to_ergs*10.2/kb/T
    an = g1/g2*np.exp(eta)*(1 + np.exp(-1.0*eta)*(W - 1.0))/W
    plt.plot(r,an,linewidth=4,color='red')

    plt.ylim(0,110)
    plt.title('2level_nlte test')
    plt.ylabel(r'level population ratio: $n_1/n_2$',size=15)
    plt.xlabel('radius/r_phot',size=15)
    plt.legend(['Sedona result','analytic result'],loc=2)

    max_err,mean_err = get_error(n1/n2,an)
    if (mean_err > 0.03): failure = 1
    if (max_err > 0.04):  failure = 1
    print 'max_error  = ',max_err
    print 'mean_error = ',mean_err

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
