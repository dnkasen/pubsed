import os
import numpy as np
import matplotlib.pyplot as plt
import h5py

def run_test(pdf="",runcommand=""):
 
    ###########################################
    # clean up old results and run the code
    ###########################################
    if (runcommand != ""): 
        os.system("rm spectrum_* plt_* integrated_quantities")
        os.system(runcommand)

     
    ###########################################
    # compare the output
    ###########################################
    plt.clf()

    h   = 6.6260755e-27    # planck's constant (ergs-s)
    c   = 2.99792458e10    # speed of light (cm/s)
    k   = 1.380658e-16     # boltzmann constant (ergs/K)
    sb  = 5.6704e-5        # stefan boltzman constant (ergs cm^-2 s^-1 K^-4) 
    pi  = 3.14159          # just pi
    plt.ion()
    
    # compare spectrum 
    plt.clf()
    fin = h5py.File('plt_00001.h5','r')
    nu = np.array(fin['nu'])
    Jnu = np.array(fin['zonedata/0/Jnu'])
    opac = np.array(fin['zonedata/0/opacity'])
    emis = np.array(fin['zonedata/0/emissivity'])
    Snu = np.nan_to_num(emis/opac)
    fin.close()
    plt.plot(nu,Jnu,'o')

    # blackbody spectrum
    T = 1e5
    f = 2.0*h*np.power(nu, 2. * np.ones_like(nu))/c**2/(np.exp(np.minimum(h*nu/k/T,35. * np.ones_like(nu))) - 1) * nu # last factor of nu added here to avoid overflow
    plt.plot(nu,f,color='red',linewidth=2)
    
    plt.legend(['sedona','analytic blackbody'])
    plt.title('one zone eq test: radiation field')
    plt.xlabel('frequency (Hz)')
    plt.ylabel('Flux')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(1e14,5e16)
    plt.ylim(1e-4,5)

    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = raw_input()


    plt.clf()

    fin = h5py.File('plt_00001.h5','r')
    b = np.array(fin['zonedata/0/Z_1/level_departure'])
    plt.plot(b,'o',color='black')
    plt.plot([0,len(b)+1],[1,1],color='black')
    plt.title('one zone eq test: NLTE departure coefficients')
    plt.xlabel('level #',size=15)
    plt.ylabel('departure coefficient',size=15)

    for i in range(len(b)):
        if (b[i] > 1.2 or b[i] < 0.8):
            plt.plot([i,i],[-1,10],color='red')

    plt.ylim(0.8,1.2)


    fin.close()
    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = raw_input()


    # this should return !=0 if failed
    return 0



if __name__=='__main__': run_test()
