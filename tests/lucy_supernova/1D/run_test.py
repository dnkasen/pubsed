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

    x,y,c = np.loadtxt('optical_spectrum_final.dat',unpack=1,skiprows=1)
    x = x/3600.0/24.0
    plt.plot(x,y,'o',color='black',markersize=8,markeredgewidth=2,markerfacecolor='none')

    # gamma-ray deposition
    tdep = []
    gdep = []
    for j in range(1,150):
        ray = 'plt_00'
        if (j< 100): ray = ray + '0' 
        if (j < 10): ray = ray + '0' 
        ray = ray + str(j) + '.dat'
        fin = open(ray,'r')
        line = fin.readline()
        tdep.append(float(line.split()[3])/3600.0/24.0)
        sum = 0
        r0 = 0 
        line = fin.readline()
        for line in fin:
            data = line.split()
            r1 = float(data[0])
            vol = 4.0*3.14159/3.0*(r1**3 - r0**3)
            sum += float(data[5])*vol
            r0 = r1
        gdep.append(sum)
        
    plt.plot(tdep,gdep,'o')
        
    # benchmark results
    x,y = np.loadtxt('../comparefiles/lucy_lc.dat',unpack=1)
    plt.plot(x,y,color='red',linewidth=2)
    x,y = np.loadtxt('../comparefiles/lucy_gr.dat',unpack=1)
    plt.plot(x,y,color='blue',linewidth=2)
    plt.ylim(1e40,0.4e44)

    plt.title
    plt.legend(['sedona LC','sedona GR','lucy LC','lucy GR'])
    plt.xlim(0,60)

    if (pdf != ''): pdf.savefig()
    else:
        plt.ion()
        pltx.show()
        j = raw_input()

    # this should return !=0 if failed
    return 0
        



if __name__=='__main__': run_test('')

