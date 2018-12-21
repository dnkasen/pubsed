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

    m_p = 1.67262158e-24   # mass of proton (g)
    k   = 1.380658e-16     # boltzmann constant (ergs/K)

    # read exact solution and mine
    z,x,rho,p,v = np.loadtxt('exact_solution.dat',unpack=True)
    r_me,rho_me,v_me,T_me = np.loadtxt('plt_00005.dat',unpack=True,usecols=[0,1,2,3])
    r_me = r_me - 1e3
    p_me = rho_me/m_p*k*T_me
    t = p*m_p/rho/k

    xlims = (0,6)
    use = (x < 5)
    plt.suptitle('sod shock tube hydro test')

    plt.subplot(2,2,1)
    plt.xlim(xlims)
    plt.plot(x,rho,color='red',linewidth=3)
    plt.plot(r_me,rho_me,'o',markerfacecolor='none')
    plt.ylabel('density')
    max_err,mean_err = get_error(rho_me,rho,x=r_me,x_comp=x,use = use)
    if (mean_err > 0.05): failure = 1

    plt.subplot(2,2,2)
    plt.xlim(xlims)
    plt.plot(x,v,color='red',linewidth=3)
    plt.plot(r_me,v_me,'o',markerfacecolor='none')
    plt.ylabel('velocity')
    max_err,mean_err = get_error(v_me,v,x=r_me,x_comp=x,use = use)
    if (mean_err > 0.05): failure = 2

    plt.subplot(2,2,3)
    plt.xlim(xlims)
    plt.plot(x,p,color='red',linewidth=3)
    plt.plot(r_me,p_me,'o',markerfacecolor='none')
    plt.ylabel('pressure')
    max_err,mean_err = get_error(p_me,p,x=r_me,x_comp=x,use = use)
    if (mean_err > 0.05): failure = 3

    plt.subplot(2,2,4)
    plt.xlim(xlims)
    plt.plot(x,t,color='red',linewidth=3)
    plt.plot(r_me,T_me,'o',markerfacecolor='none')
    plt.ylabel('temperature')
    max_err,mean_err = get_error(T_me,t,x=r_me,x_comp=x,use = use)
    if (mean_err > 0.05): failure = 4

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
