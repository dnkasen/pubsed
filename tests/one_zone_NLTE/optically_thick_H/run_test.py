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
    	os.system("rm spectrum_* plt_* integrated_quantities")
    	os.system(runcommand)

    fail_flag = 0
    ###########################################
    # compare the output
    ###########################################
    plt.clf()


    h   = 6.6260755e-27    # planck's constant (ergs-s)
    c   = 2.99792458e10    # speed of light (cm/s)
    k   = 1.380658e-16     # boltzmann constant (ergs/K)
    sb  = 5.6704e-5        # stefan boltzman constant (ergs cm^-2 s^-1 K^-4) 
    pi  = 3.14159          # just pi
    #    plt.ion()

    T   = 5.e4
    data = np.loadtxt('plt_00003.dat',skiprows=2)

    #------------------------------------------
    #compare output spectrum
    #------------------------------------------

    plt.clf()
    # compare spectrum 

    fin = h5py.File('plt_00003.h5','r')
    nu = np.array(fin['nu'])
    Jnu = np.array(fin['zonedata/0/Jnu'])
    fin.close()
    plt.plot(nu,Jnu,'o',color='black')

    # blackbody spectrum
    f = 2.0*h*np.power(nu,2.0 * np.ones_like(nu))/pow(c,2.)/(np.exp(np.minimum(h*nu/k/T,60. * np.ones_like(nu))) - 1.)*nu # written like this to avoid overflow
    plt.plot(nu,f,color='red',linewidth=2)

    plt.legend(['sedona','analytic blackbody'])
    plt.title('one zone NLTE test for optically thick H: radiation field')
    plt.xlabel('frequency (Hz)')
    plt.ylabel('Flux')
    plt.yscale('log')
    #py.xscale('log')
    plt.xlim(0,2.5e16)
    plt.ylim(1.e-8,1.e-1)


    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = raw_input()



    #------------------------------------------
    #check LTE departure coefficients
    #------------------------------------------
    plt.clf()

    
    plt.title('one zone NLTE test for optically thick H: LTE departure coefficients')
    plt.xlabel('level #',size=15)
    plt.ylabel('departure coefficient',size=15)

    plt.plot([0,35],[1,1],color='red')

    fin = h5py.File('plt_00003.h5','r')
    b = np.array(fin['zonedata/0/Z_1/level_departure'])
    fin.close()
    plt.plot(b,'o',color='black')
    

    for i in range(len(b)):
        if (b[i] > 1.1 or b[i] < 0.9):
            plt.plot([i,i],[-1,10],color='red')
            fail_flag = 1

    #    plt.ylim(0.5,1.5)
    plt.xlim(0,35)



    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = raw_input()

    # this should return !=0 if failed
    return fail_flag



if __name__=='__main__': run_test('')

