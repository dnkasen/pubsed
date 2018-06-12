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

    # physics constants
    h   = 6.6260755e-27    # planck's constant (ergs-s)
    c   = 2.99792458e10    # speed of light (cm/s)
    k   = 1.380658e-16     # boltzmann constant (ergs/K)
    sb  = 5.6704e-5        # stefan boltzman constant (ergs cm^-2 s^-1 K^-4) 
    pi  = 3.14159          # just pi

    # parameters of the model
    testname = "2D radiating sphere"
    T   = 1.0e4
    L  = 1e43
    r0 = 0.5e15

    # open up and read output file
    fin = h5py.File('plt_00001.h5','r')
    dr   = np.array(fin['dr'])
    z    = np.array(fin['z'])
    p    = np.array(fin['r'])
    Trad = np.array(fin['T_rad'])


    #------------------------------------------
    # plot temperature image
    #------------------------------------------
    plt.matshow(Trad)    
    plt.title(testname + ': radiation field')
    
    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = raw_input()
    plt.clf()


    
    #------------------------------------------
    # compare 1D slice of temperature
    #------------------------------------------
    for j in (25,50,100):
        plt.plot(p,Trad[:,j],'o')

        # dillution factor W
        rr = (p*p + z[j]*z[j])**0.5
        w    = 0.5*(1 - (1 - r0**2/rr**2)**0.5)
        TW = (L*w/(4.0*pi*r0**2)/sb)**0.25
        plt.plot(p,TW,color='k',linewidth=3)
    
    plt.legend(['sedona Trad','analytic solution'])
    plt.title(testname + ': radiation field')
    plt.xlabel('radius (cm)')
    plt.ylabel('radiation temperature (aT^4 = erad)')

    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = raw_input()
    fin.close()

    #------------------------------------------
    # compare output spectrum
    #------------------------------------------
        
    plt.clf()
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
        j = raw_input()


    # this should return !=0 if failed
    return 0
        



if __name__=='__main__': run_test('')

