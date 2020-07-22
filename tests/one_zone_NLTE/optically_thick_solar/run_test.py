import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys




def run_test(pdf="",runcommand=""):

    fail_flag = 0

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
    plt.title('one zone NLTE test for optically thick solar elements: radiation field')
    plt.xlabel('frequency (Hz)')
    plt.ylabel('Flux')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(2.5e13,2.5e17)
    plt.ylim(1.e-10,1.e-1)

    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = raw_input()

    #------------------------------------------
    #check LTE departure coefficients
    #------------------------------------------
    plt.clf()

    elements = []
    element_ion_stage_indices = []

    for groupname in fin['/zonedata/0'].keys():
        if (groupname[0] == 'Z'):
            elements.append(str(groupname))
            ion_stage_indices = []
            max_ion_stages = len(fin['/zonedata/0/' + str(groupname) + '/ion_fraction'] )
            sample_ion_fractions = np.array((fin['/zonedata/0/' + str(groupname) + '/ion_fraction'] ))
            for ion_stage in range(max_ion_stages):
                if (sample_ion_fractions[ion_stage] > -1):
                    ion_stage_indices.append(ion_stage)
            element_ion_stage_indices.append(ion_stage_indices)

    num_elements = len(elements)
        
    Z_elem = []
    for elem in elements:
        Z_elem.append(int(elem.split('_')[1]))


    nrow = num_elements
    ncol = 1
    fig, axs = plt.subplots(nrows=nrow, ncols=ncol)

    lte_fraction_threshhold_plot = 1.e-11
    lte_fraction_threshhold_test = 1.e-9
    fig.suptitle("Departure coefficients, optically thick solar\nNot plotting levels with LTE fraction < %.1e\nTest success depends only on levels with LTE fraction > %.1e\n\n" % (lte_fraction_threshhold_plot,lte_fraction_threshhold_test), fontsize=10)


    atomic_database =h5py.File('./atom_data.hdf5', 'r')
    element_index = 0
    for ax in axs.reshape(-1):
        this_Z = str(Z_elem[element_index])
        num_stages = len(element_ion_stage_indices[element_index])
        ion_ground_states = np.array(atomic_database[str(this_Z) + '/ion_ground'])
        level_bs = np.array(fin['zonedata/0/Z_' + this_Z + '/level_departure'])
        level_ns = np.array(fin['zonedata/0/Z_' + this_Z + '/level_fraction'])
        level_n_ltes = level_ns / level_bs 
        num_levels_element = len(level_bs)
        for level_index in range(num_levels_element):
            if (level_n_ltes[level_index] > lte_fraction_threshhold_plot):
                if (level_bs[level_index] <= 1.5 and level_bs[level_index]  >= 0.5):
                    ax.plot(level_index, level_bs[level_index],'o',color='green')
                elif (level_bs[level_index] >= 1.5):
                    if (level_n_ltes[level_index] >= lte_fraction_threshhold_test):
                        fail_flag = 1
                    ax.plot(level_index,1.5,'o',c='r')
                else:
                    if (level_n_ltes[level_index] >= lte_fraction_threshhold_test):
                        fail_flag = 1
                    ax.plot(level_index,0.5,'o',c='r')
            else:
                ax.axvline(level_index,0,2,linewidth=2,c=plt.cm.binary(50))
        ax.plot([0,num_levels_element],[1,1],color='black')
        for ground_state in ion_ground_states:
            ax.axvline(ground_state,0,2,linewidth=2,c='k')
        ax.set_ylim(0.5,1.5)
        ax.set_title("Z = " + this_Z,fontsize=8)
        element_index += 1
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.82)

    fin.close()
    atomic_database.close()


    if (pdf != ''): pdf.savefig()
    else:
        plt.show()
        j = raw_input()

    # this should return !=0 if failed
    plt.clf()
    return fail_flag


if __name__=='__main__': run_test('')
