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
     
    ###########################################
    # compare the output
    ###########################################
    plt.clf()
    plt.ion()


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

    plt.plot(nu,Jnu,'o',color='black')

    # blackbody spectrum
    f = 2.0*h*np.power(nu,2.0 * np.ones_like(nu))/pow(c,2.)/(np.exp(np.minimum(h*nu/k/T,60. * np.ones_like(nu))) - 1.)*nu # written like this to avoid overflow
    plt.plot(nu,f,color='red',linewidth=2)

    plt.legend(['sedona','analytic blackbody'])
    plt.title('one zone NLTE test for optically solar elements: radiation field')
    plt.xlabel('frequency (Hz)')
    plt.ylabel('Flux')
    plt.yscale('log')
    plt.xscale('log')
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

    with open("onezone.mod", 'r') as modelfile:
        for i, line in enumerate(modelfile):
            if (i == 2):
                elem_id_string = line

    Z_elem = []
    for elem_id in elem_id_string.split():
        Z_elem.append(elem_id.split('.')[0])

    num_elems = len(Z_elem)
    counter = 0

    atomic_database = h5py.File(os.environ['SEDONA_HOME'] + str('/data/atom_data.hdf5'), 'r')

    nrow = num_elems
    ncol = 1
    fig, axs = plt.subplots(nrows=nrow, ncols=ncol)

    fig.suptitle("Departure coefficients, optically thin solar\n\n", fontsize=10)

    counter = 0
    for ax in axs.reshape(-1):
        this_Z = str(Z_elem[counter]) 
        b = np.array(fin['zonedata/0/Z_' + this_Z + '/level_departure'])
        ax.plot([-1 * len(b)/10,len(b) + len(b)/10],[1,1],color='black')
        ax.plot(b,'o',color='green')
        ion_ground_states = np.array(atomic_database[this_Z + '/ion_ground'])
        for ground_state in ion_ground_states:
            ax.axvline(ground_state,0,2)
        ax.set_ylim(0.5,1.5)
        ax.set_title("Z = " + this_Z,fontsize=8)
        counter += 1

    plt.tight_layout()
    plt.subplots_adjust(top=0.9)

    fin.close()
    atomic_database.close()


    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = raw_input()

    # this should return !=0 if failed
    plt.clf()
    return 0


if __name__=='__main__': run_test('')

