import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys
import bisect


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
    # compare the bolometric light curve
    #-------------------------------------------
    failure = 0
    plt.clf()

    # read the spectral data
    fin = h5py.File('spectrum_final.h5','r')
    nu    = np.array(fin['nu'])
    times = np.array(fin['time'])
    Lnu   = np.array(fin['Lnu'])
    lam = 3e10/nu*1e8
    days = times/3600/24.0
    fin.close()

    # plot bolometric light curve
    nt  = len(times)
    nnu = len(nu)
    bol = np.zeros(nt)
    for it in range(nt):
        bol[it] = np.trapz(Lnu[it,:],x=nu)
    plt.plot(days,bol,color='k',linewidth=3)

    # read the comparison data
    fin = h5py.File('comparison_files/lightcurve_output.h5','r')
    nu_c    = np.array(fin['nu'])
    times_c = np.array(fin['time'])
    Lnu_c   = np.array(fin['Lnu'])
    lam_c   = 3e10/nu*1e8
    days_c  = times/3600/24.0
    fin.close()

    # plot bolometric light curve
    nt_c  = len(times_c)
    nnu_c = len(nu_c)
    bol_c = np.zeros(nt_c)
    for it in range(nt_c):
        bol_c[it] = np.trapz(Lnu_c[it,:],x=nu_c)

    plt.plot(days_c,bol_c,'o',color='r',linewidth=3)

    plt.xlim(0,60)
    plt.ylim(1e41,3e43)
    plt.xlabel('days since explosion')
    plt.ylabel('bolometric luminosity')
    plt.yscale('log')
    plt.legend(['sedona output','reference'])


    use = (days > 2)
    max_err,mean_err = get_error(bol,bol_c,use=use)
    if (mean_err > 0.1): failure = 1
    if (max_err > 0.2): failure = 1

    plt.title('1D toy Ia supernova - bolometric light curve')

    # add plot to pdf file (or show on screen)
    if (pdf != ''): pdf.savefig()
    else:
        plt.ion()
        plt.show()
        j = get_input('Press any key to continue >')



    #-------------------------------------------
    # compare the spectra
    #-------------------------------------------
    plt.clf()

    plt.title('1D toy Ia supernova - spectra')

    tplot = [5,10,20,40]
    cnt = 1
    for tc in tplot:
        plt.subplot(2,2,cnt)

        # plot results
        indt = bisect.bisect(days,float(tc))
        Llam = Lnu[indt,:]*nu**2/2.99e10/1e8
        plt.plot(lam/1e4,Llam/1e39,color='k')
        # plot comparison
        indt = bisect.bisect(days_c,float(tc))
        Llam_c = Lnu_c[indt,:]*nu_c**2/2.99e10/1e8
        plt.plot(lam_c/1e4,Llam_c/1e39,color='r')

        plt.xlim(0.1,1.2)
        plt.annotate('day ' + str(tc),(0.8,max(Llam_c/1e39)*0.8))
        if (cnt == 1):
            plt.legend(['sedona output','reference'],loc=4)

        # check errors
        use = (Llam_c > max(Llam_c)*0.1)
        max_err,mean_err = get_error(Llam,Llam_c,use=use)
        if (mean_err > 0.1): failure = 2

        #    plt.xlim(1000,10000)
        #    plt.xlabel('wavelength (angstroms)')
        #    plt.ylabel('specific luminoisty (ergs/s/angstrom)')
        #    plt.title(testname)
        cnt += 1

    # add plot to pdf file (or show on screen)
    if (pdf != ''): pdf.savefig()
    else:
        plt.ion()
        plt.show()
        j = get_input('Press any key to continue >')


        #    Llam = Lnu*nu/lam
        #    plt.plot(lam,Llam,color='k',lw=2)

            # plot the comparison spectrum
        #    nu,Lnu = np.loadtxt('reference_spectrum.dat',unpack=True,skiprows=1,usecols=[0,1])
        #    lam = 3e10/nu*1e8
        #    Llam_ref = Lnu*nu/lam
        #    plt.plot(lam,Llam_ref,color='r',lw=2)

        #    plt.xlim(1000,10000)
        #    plt.xlabel('wavelength (angstroms)')
        #    plt.ylabel('specific luminoisty (ergs/s/angstrom)')
        #    plt.legend(['sedona output','reference'])
        #    plt.title(testname)


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
