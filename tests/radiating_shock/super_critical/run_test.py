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

    testname = "radiating_subcritical_shock"

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

    for i in [9]:
        name = 'plt_0000' + str(i) + '.h5'
        if (i > 9): name = 'plt_000' + str(i) + '.h5'
        print (name)
        fin = h5py.File(name,'r')
        time = (np.array(fin['time']))[0]
        r  = np.array(fin['r'])
        T_gas = np.array(fin['T_gas'])
        T_rad = np.array(fin['T_rad'])
        print (time/1e4)
        plt.plot(r,T_gas,color='k')
        plt.plot(r,T_rad,'--',color='r')

    # comparisons -- output every 5e2 seconds
#    r,T_gas,T_rad = np.loadtxt('comparison/ray_00004',usecols=[0,3,4],unpack=True)
 #   plt.plot(r-2.4e12 + 8e11,T_gas,color='b')
#   plt.plot(r-2.4e12 + 8e11,T_rad,'--',color='b')
  #  r,T_gas,T_rad = np.loadtxt('comparison/ray_00008',usecols=[0,3,4],unpack=True)
   # plt.plot(r-2.4e12 + 8e11,T_gas,color='b')
#    plt.plot(r-2.4e12 + 8e11,T_rad,'--',color='b')
#    
    r,T_gas,T_rad = np.loadtxt('comparison_files/ray_00200_zeus',usecols=[0,4,5],unpack=True)
    T_rad = T_rad*11614.5300726
    T_gas = T_gas*11614.5300726
    r = 2e6*4.0e3 + r + 8e11
    plt.plot(r,T_gas)
    plt.plot(r,T_rad)



    plt.xlabel('radius')
    plt.ylabel('temperature')
    plt.xlim(8.0e11,8.4e11)
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
