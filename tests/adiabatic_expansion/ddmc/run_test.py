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

    t0    = 20*3600*24
    vedge = 0.5e9
    temp0 = 1e4
    npts = 50
    time = np.zeros(npts)
    Tmod = np.zeros(npts)

    colors = ['k','b','g','r','y','pink']
    cnt = 0
    for j in range(0,npts):

        name = 'plt_000' + str(j) + '.h5'
        if (j < 10):
            name = 'plt_0000' + str(j) + '.h5'
        fin = h5py.File(name,'r')
        time[j] = (np.array(fin['time']))[0]
        r = np.array(fin['r'])
        T_rad = np.array(fin['T_rad'])
        redge = vedge*time[j]

        Tmod[j] = np.mean(T_rad[(r > redge*0.1)*(r < redge*0.7)])

        if (j % 10 == 0):
            Tan = r*0 + temp0*(t0/time[j])
            Tan[r > redge] = 0
            plt.plot(r,Tan,color=colors[cnt])
            plt.plot(r,T_rad,'o',color=colors[cnt])
            cnt += 1

        plt.legend(['sedona DDMC','analytic'])
        plt.title('adiabatic expansion ddmc - temperature profiles')
        plt.xlabel('radius (cm)')
        plt.ylabel('radiation temperature')


    # add plot to pdf file (or show on screen)
    if (pdf != ''): pdf.savefig()
    else:
        plt.ion()
        plt.show()
        j = get_input('Press any key to continue >')

    plt.clf()
    plt.plot(time/3600/24.,Tmod,'o')
    Tan = temp0*(t0/time)
    plt.plot(time/3600/24.,Tan)
    plt.legend(['sedona DDMC','analytic'])
    plt.title('adiabatic expansion ddmc - mean temperature evolution')
    plt.xlabel('time (days)')
    plt.ylabel('mean temperature')


    max_err,mean_err = get_error(Tmod,Tan)
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
