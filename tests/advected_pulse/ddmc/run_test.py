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

    testname = "advected pulse"

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
	# read initial energy
    r,v,d,tr,comp = np.loadtxt('../models/constant.mod',unpack=True,skiprows=3)
    vel = v[0]
    rho = d[0]
    n = len(r)
    r0 = 1e5

    # total internal energy in pulse
    esum = 0.0
    for i in range(1,n-1):
        esum += tr[i]**4*(r[i]-r[i-1])


    colors = ['k','r','b','g']
    cnt = 0
    for i in np.arange(5,21,5):

        # read sedona data
        fname = 'plt_000' + str(i) + '.h5'
        if (i < 10): fname = 'plt_0000' + str(i) + '.h5'
        fin = h5py.File(fname,'r')
        r = np.array(fin['r'])
        T = np.array(fin['T_rad'])
        time = (np.array(fin['time']))[0]
        fin.close()

        # plot up
        plt.plot(r-r0,T,'o',color=colors[cnt],markersize=5,markeredgewidth=2,markerfacecolor='none')

        # analytic solution
        D = 3e10/rho*4.0/3.0
        T0 = (esum/(D*time*3.141519)**0.5)**0.25
        rc = r0 + time*vel
        Tan = T0*np.exp(-1.0*(r-rc)**2/(4.0*D*time))
        plt.plot(r - r0,Tan,color=colors[cnt],lw=3)
        cnt += 1

        # determine error
        use = (T > T*0.1)
        max_err,mean_err = get_error(Tan,T)
        #print('mean error = ' + str(mean_err))
        #print('max error  = ' + str(max_err))

        if (mean_err > 0.1): failure = 1
        if (np.isnan(mean_err)): failure = 1

    plt.title('Advected Pulse - Discrete Diffusion MC')
    plt.legend(['Sedona DDMC','analytic diffusion'])
    plt.xlabel('radius (cm)',size=15)
    plt.ylabel('temperature (K)',size=15)


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
