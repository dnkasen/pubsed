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
    failure = 0
    plt.clf()

    h   = 6.6260755e-27    # planck's constant (ergs-s)
    c   = 2.99792458e10    # speed of light (cm/s)
    k   = 1.380658e-16     # boltzmann constant (ergs/K)
    sb  = 5.6704e-5        # stefan boltzman constant (ergs cm^-2 s^-1 K^-4) 
    pi  = 3.14159          # just pi
    plt.ion()

    T   = 1e4
    L  = 1e43
    r0 = 0.5e15

    data = np.loadtxt('plt_00001.dat',skiprows=2)
    r    = data[:,0]
    r    = r - 0.5*(r[1] - r[0])
    tgas = data[:,3]
    trad = data[:,4]

    # dillution factor W
    rr = np.arange(r0,max(r),0.05*r0)
    w    = 0.5*(1 - (1 - r0**2/rr**2)**0.5)
    TW = (L*w/(4.0*pi*r0**2)/sb)**0.25
    w    = 0.5*(1 - (1 - r0**2/r**2)**0.5)
    tan = (L*w/(4.0*pi*r0**2)/sb)**0.25

    # do numerical comparisons
    max_err,mean_err = get_error(tgas,tan,use = (r> 5e14))
    if (max_err > 0.01): failure = 1
    max_err,mean_err = get_error(trad,tan,use = (r> 5e14))
    if (max_err > 0.01): failure = 2

    plt.plot(r,trad,'o',color='black')
    plt.plot(r,tgas,'--',color='blue')
    plt.plot(rr,TW,color='red',linewidth=2)
    plt.legend(['sedona Trad','sedona Tgas','analytic solution'])
    plt.title('spherical lightbulb test: radiation field')
    plt.xlabel('radius (cm)')
    plt.ylabel('radiation temperature (aT^4 = erad)')


    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = raw_input()
        

    #------------------------------------------
    #compare output spectrum
    #------------------------------------------

    plt.clf()
    data = np.loadtxt('spectrum_1.dat')
    nu = data[:,0]
    y  = data[:,1]
    plt.plot(nu,y,'o',color='black')
    # blackbody spectrum
    f = 2.0*h*nu**3/c**2/(np.exp(h*nu/k/T) - 1)
    f = f/(sb*T**4/pi)*L
    plt.plot(nu,f,color='red',linewidth=2)

    # do numerical comparisons
    max_err,mean_err = get_error(y,f)
    if (mean_err > 0.1): failure = 3

    plt.legend(['sedona','analytic blackbody'])
    plt.title('spherical lightbulb test: output spectrum')
    plt.xlabel('frequency (Hz)')
    plt.ylabel('Flux')
    plt.yscale('log')
    plt.xlim(0,3e15)
    plt.ylim(1e25,1e29)

    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = raw_input()

    #------------------------------------------
    # compare spectrum at zone 50
    #------------------------------------------

    plt.clf()
    fin  = h5py.File('plt_00001.h5','r')
    nu   = np.array(fin['nu'])
    Jnu  = np.array(fin['zonedata/50/Jnu'])
    rz   = np.array(fin['r'])
    fin.close()
    plt.plot(nu,Jnu,'o')

    # blackbody spectrum
    T = 1e4
    L  = 1e43
    f = 2.0*h*nu**2.0/c**2/(np.exp(h*nu/k/T) - 1)*nu
    f = L*f/(sb*T**4*4.0*pi*r0**2)
    W = 0.5*(1 - (1 - (r0/rz[50])**2)**0.5)
    f = W*f
    plt.plot(nu,f,color='red',linewidth=2)
    
    # do numerical comparisons
    max_err,mean_err = get_error(Jnu,f)
    if (mean_err > 0.1): failure = 4

    plt.legend(['sedona','analytic blackbody'])
    plt.title('spherical lightbulb test: Jnu at zone=50')
    plt.xlabel('frequency (Hz)')
    plt.ylabel('Flux')
    plt.yscale('log')
    plt.xlim(0,3e15)
    plt.ylim(5e-8,5e-4)

    
    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = raw_input()

    # this should return !=0 if failed
    return failure
        
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

if __name__=='__main__': run_test('')



