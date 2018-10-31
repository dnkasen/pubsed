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
    plt.ion()
    failure = 0
    testname = "2D spherical lightbulb"

    # physics constants
    h   = 6.6260755e-27    # planck's constant (ergs-s)
    c   = 2.99792458e10    # speed of light (cm/s)
    k   = 1.380658e-16     # boltzmann constant (ergs/K)
    sb  = 5.6704e-5        # stefan boltzman constant (ergs cm^-2 s^-1 K^-4)
    pi  = 3.14159          # just pi

    # parameters of the model
    T   = 1.0e4
    L  = 1e43
    r0 = 0.5e15

    # open up and read output file
    fin = h5py.File('plt_00001.h5','r')
    dr   = np.array(fin['dr'])
    z    = np.array(fin['z'])
    p    = np.array(fin['r'])
    Trad = np.array(fin['T_rad'])

    # shift to cell centered
    p = p + dr[0]/2.0
    z = z + dr[1]/2.0

    #------------------------------------------
    # plot temperature image
    #------------------------------------------
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.matshow(Trad)

    ax.set_title(testname + ': radiation field')

    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = get_input('Press any key to continue >')
    plt.clf()


    #------------------------------------------
    # compare 1D slices of temperature
    #------------------------------------------
    for j in ([25,50,100]):

        thisT = Trad[:,j]
        plt.plot(p,thisT,'o')

        # Get analytic temeprature from dillution factor W
        rr = (p*p + z[j]*z[j])**0.5
        # Make rr < r0 a zero temperature (i.e., rr --> infinty)
        rr[rr < r0] = r0*1e10
        w    = 0.5*(1 - (1 - r0**2/rr**2)**0.5)
        TW = (L*w/(4.0*pi*r0**2)/sb)**0.25
        plt.plot(p,TW,color='k',linewidth=3)

        use = []
        if (j == 100): use = (p > 5.1e14)
        max_err,mean_err = get_error(thisT,TW,use=use)
        if (max_err > 0.05 or mean_err > 0.01):
            failure = 1

    plt.legend(['sedona Trad','analytic solution'])
    plt.title(testname + ': radiation field')
    plt.xlabel('radius (cm)')
    plt.ylabel('radiation temperature (aT^4 = erad)')

    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = get_input('Press any key to continue >')
    fin.close()
    plt.clf()

    #------------------------------------------
    # compare output spectrum
    #------------------------------------------

    fin = h5py.File('spectrum_1.h5','r')
    nu  = np.array(fin['nu'])
    Lnu = np.array(fin['Lnu'])
    nmu = (Lnu.shape)[2]
    for i in range(nmu):
        plt.plot(nu,Lnu[0,:,i],'o') #,color=cmap[i])

    # blackbody spectrum
    f = 2.0*nu*h*nu**2.0/c**2/(np.exp(h*nu/k/T) - 1)
    f = f/(sb*T**4/pi)*L
    plt.plot(nu,f,color='black',linewidth=5)
    plt.plot(nu,f,color='white',linewidth=3)

    # calculate error
    for i in range(nmu):
        max_err,mean_err = get_error(Lnu[0,:,i],f)
        if (mean_err > 0.1): failure = 2

    plt.legend(['sedona','analytic blackbody'])
    plt.title(testname + ': output spectrum')
    plt.xlabel('frequency (Hz)')
    plt.ylabel('Flux')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(2e13,4e15)
    plt.ylim(1e24,3e28)

    if (pdf != ''): pdf.savefig()
    else:
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
