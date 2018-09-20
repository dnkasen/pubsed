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
    	os.system("rm spectrum_* plt_* integrated_quantities.dat")
    	os.system(runcommand)
     
    #-------------------------------------------
    # compare the output
    #-------------------------------------------
    failure = 0
    plt.clf()

    # write code here to read in output data and 
    # plot it compared to standard reference

#    for example...
#    plt.plot(x,y)
#    plt.legend(['sedona','analytic'])
#    plt.title(testname)
#    plt.xlabel('x')
#    plt.ylabel('y')
#    plt.yscale('log')
#    plt.xlim(0,1)
#    plt.ylim(0,1)


    # write code here to do a numerical comparison
    # of output data to standard reference and
    # determine success or failure
    # set failure to some number != 0 if something
    # is wrong

#    max_err,mean_err = get_error(y_sedona,y_reference)
#    if (mean_err > 0.1): failure = 1


    # add plot to pdf file (or show on screen)
    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = raw_input()

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

if __name__=='__main__': run_test('')



