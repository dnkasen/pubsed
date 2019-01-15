import numpy as np
import h5py

####################################
nr   = 1
rho0    = 1.e-7
Tgas0     = 5.5e8
erad = 1.e12 # erg per cm^3
rmin = 0.
rmax  = 1.e10 # arbitrary
vmax   = 0.
t0   = 0.
Z  = [1]
A  = [1.]
mass_fracs = [1.]
###################################

##################################
# Make sedona 1D spherical model
#################################

comp = np.zeros((nr,len(Z)))
for i in range(nr):
    for k in range(len(Z)):
        comp[i,k] = mass_fracs[k]
	
# in general would create arrays, loop through zones. But here just using one zone
fout = h5py.File('onezone.h5','w')
fout.create_dataset('time',data=[t0],dtype='d')
fout.create_dataset('r_min',data=[rmin],dtype='d')
fout.create_dataset('r_out',data=[rmax],dtype='d')
fout.create_dataset('Z',data=Z,dtype='i')
fout.create_dataset('A',data=A,dtype='i')
fout.create_dataset('rho',data=[rho0],dtype='d')
fout.create_dataset('temp',data=[Tgas0],dtype='d')
fout.create_dataset('v',data=[vmax],dtype='d')
fout.create_dataset('erad',data=[erad],dtype='d')
fout.create_dataset('comp',data=comp,dtype='d')


