import matplotlib
from matplotlib.pyplot import *
import numpy as np
from math import exp, sqrt, log, log10
import h5py


###########################################
# compare the output
###########################################



h   = 6.6260755e-27    # planck's constant (ergs-s)
c   = 2.99792458e10    # speed of light (cm/s)
k   = 1.380658e-16     # boltzmann constant (ergs/K)
sb  = 5.6704e-5        # stefan boltzman constant (ergs cm^-2 s^-1 K^-4) 
pi  = 3.14159          # just pi
arad = 4. * sb / c


fin = h5py.File('plt_00001.h5','r')
nu = np.array(fin['nu'])

zone_information = np.loadtxt("plt_00000.dat")
initial_x = 0.01
Te = zone_information[3]

x_frequencies = np.fromiter( (h * nu[i] / (k * Te) for i in range(len(nu))),np.float)

initial_energy_density = arad * pow(zone_information[4],4.)

initial_frequency = initial_x * k * Te / h

initial_photon_number_density = initial_energy_density / (h * initial_frequency)
N_gamma = initial_photon_number_density * pow(c,3.)/(8. * np.pi)
chemical_potential = -1. * log(1./2. * pow(h/(k * Te),3.) * N_gamma  )
print "chemical potential is " + str(chemical_potential)

equilibrium_occupation_distribution = np.fromiter( (1./(exp(chemical_potential + x_frequencies[index])) for index in range(len(x_frequencies))),np.float)

equilibrium_Jnu = np.fromiter( (2. * equilibrium_occupation_distribution[index] * h * pow(k * Te * x_frequencies[index]/h,3.) /pow(c,2.) for index in range(len(x_frequencies))),np.float)


BE_occupation_distribution = np.fromiter( (1./(exp(chemical_potential + x_frequencies[index]) - 1.) for index in range(len(x_frequencies))),np.float)
BE_Jnu = np.fromiter( (2. * BE_occupation_distribution[index] * h * pow(k * Te * x_frequencies[index]/h,3.) /pow(c,2.) for index in range(len(x_frequencies))),np.float)


ntime = 30
for i in range(0,ntime,2):

    if i < 10:
        time_string = '0000' + str(i)
    else:
        time_string = '000' + str(i)

    #handle first time as special case
    if i == 0:
        time_string = '00001'
    
    fin = h5py.File('plt_' + time_string + '.h5','r')
    nu = np.array(fin['nu'])
    Jnu = np.array(fin['zonedata/0/Jnu'])
    rz = np.array(fin['r'])
    fin.close()
    loglog(x_frequencies, Jnu,'-',markersize=3,linewidth=2,c=cm.RdYlGn(255 - i * 255/ntime))


loglog(x_frequencies,equilibrium_Jnu,'--',linewidth=3,c='k')
#loglog(x_frequencies,BE_Jnu,':',linewidth=3,c='m')

xlabel(r"$h \nu / k_B T_e$",fontsize=24)
xlim(3.e-4,30.)
ylim(1.e-4,2.e-1)

ylabel(r"$J_{\nu} \,\, {\rm (erg \,\,cm}^{-2}\,\,{\rm s}^{-1} \, \, {\rm hz}^{-1}{\rm)}$",fontsize=24)

leg = legend(fontsize=14,framealpha=1.)

ax = subplot(111)

ax.tick_params(axis='both', which='major', labelsize=18)

show()
close()

cohr = []
ntime = 30
for i in range(0,ntime,1):

    if i < 10:
        time_string = '0000' + str(i)
    else:
        time_string = '000' + str(i)

    #handle first time as special case
    if i == 0:
        time_string = '00001'
    
    fin = h5py.File('plt_' + time_string + '.h5','r')
    this_cohr = np.array(fin['compton_heating_rate'])
    fin.close()
    cohr.append(this_cohr)

plot(cohr)
show()
close()

