import numpy as np
import h5py
import matplotlib.pyplot as plt

####################################
m_sun  = 1.99e33
pi     = 3.14159
mass   = 1.4*1.99e33
vmin   = 0
vmax   = 1.0e9
texp   = 1.0*(3600.0*24.0)
Z  = [14,26,27,28]
A  = [28,56,56,56]
T0     = 1.0e4*(20)
###################################



rmin0   = vmin*texp
rmax0   = vmax*texp
rho0    = mass/((4.0*pi/3.0)*(rmax0**3))
n_elems = len(Z)



##################################
# Make sedona 3D hdf5 model
#################################
nr = 100
r_out = np.linspace(rmin0, rmax0, nr+1, endpoint=True)[1:]
vr = r_out/texp
nr_array = np.arange(nr)

ntheta = 10
costheta_out = np.linspace(1, -1, ntheta+1, endpoint=True)[1:]
theta_out = np.arccos(costheta_out)
# theta_out = np.linspace(0, np.pi, ntheta+1, endpoint=True)[1:]
vtheta = np.zeros(ntheta)
ntheta_array = np.arange(ntheta)

nphi = 10
phi_out = np.linspace(0, 2*np.pi, nphi+1, endpoint=True)[1:]
vphi = np.zeros(nphi)
nphi_array = np.arange(nphi)

rho  = np.zeros((nr,ntheta,nphi))
temp = np.zeros((nr,ntheta,nphi))
comp = np.zeros((nr,ntheta,nphi,len(Z)))
erad = np.zeros((nr,ntheta,nphi))
velr = np.zeros((nr,ntheta,nphi))
veltheta = np.zeros((nr,ntheta,nphi))
velphi = np.zeros((nr,ntheta,nphi))

for i in range(nr):
	for j in range(ntheta):
		for k in range(nphi):

			rho[i,j,k] = rho0
			temp[i,j,k] = T0

			velr[i,j,k] = vr[i];
			veltheta[i,j,k] = vtheta[j];
			velphi[i,j,k] = vphi[k];

			# get composition
			m_enc = (4.0*pi/3.0)*(r_out[i]**3.0)*rho0/m_sun	
			if   (m_enc < 0.50): ni_frac = 1.0
			elif (m_enc < 0.75): ni_frac = (0.75 - m_enc)/0.25
			else: ni_frac = 0
			comp[i][j][k][0] = 1 - ni_frac
			comp[i][j][k][1] = 0.0
			comp[i][j][k][2] = 0.0
			comp[i][j][k][3] = ni_frac;

fout = h5py.File('lucy_3D-testing_spherical_grid.h5','w')
fout.create_dataset('time',data=[texp],dtype='d')
fout.create_dataset('Z',data=Z,dtype='i')
fout.create_dataset('A',data=A,dtype='i')
fout.create_dataset('rho',data=rho,dtype='d')
fout.create_dataset('temp',data=temp,dtype='d')
fout.create_dataset('vr',data=velr,dtype='d')
fout.create_dataset('vtheta',data=veltheta,dtype='d')
fout.create_dataset('vphi',data=velphi,dtype='d')
fout.create_dataset('erad',data=erad,dtype='d')
fout.create_dataset('comp',data=comp,dtype='d')
fout.create_dataset('r_out',data=r_out,dtype='d')
fout.create_dataset('theta_out',data=theta_out,dtype='d')
fout.create_dataset('phi_out',data=phi_out,dtype='d')
fout.create_dataset('rmin',data=[rmin0],dtype='d')