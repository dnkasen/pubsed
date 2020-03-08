import numpy as np
import h5py
import matplotlib.pyplot as plt

####################################
m_sun  = 1.99e33
pi     = 3.14159
mass   = 1.4*1.99e33
vmax   = 1.0e9
texp   = 1.0*(3600.0*24.0)
Z  = [14,26,27,28]
A  = [28,56,56,56]
T0     = 1.0e4*(20)
###################################



rmax    = vmax*texp
rho0    = mass/(4.0*pi/3.0*rmax**3)
n_elems = len(Z)



##################################
# Make sedona 3D hdf5 model
#################################
nx_1 = 15
vx_1 = np.linspace(-vmax, -vmax/4., nx_1+1, endpoint=True)[1:]
nx_2 = 30
vx_2 = np.linspace(-vmax/4., vmax/4., nx_2+1, endpoint=True)[1:]
nx_3 = 15
vx_3 = np.linspace(vmax/4., vmax, nx_3+1, endpoint=True)[1:]
nx = nx_1 + nx_2 + nx_3
vx = np.r_[vx_1, vx_2, vx_3]
nx_array = np.arange(nx)

ny_1 = 15
vy_1 = np.linspace(-vmax, -vmax/4., ny_1+1, endpoint=True)[1:]
ny_2 = 30
vy_2 = np.linspace(-vmax/4., vmax/4., ny_2+1, endpoint=True)[1:]
ny_3 = 15
vy_3 = np.linspace(vmax/4., vmax, ny_3+1, endpoint=True)[1:]
ny = ny_1 + ny_2 + ny_3
vy = np.r_[vy_1, vy_2, vy_3]
ny_array = np.arange(ny)

nz_1 = 15
vz_1 = np.linspace(-vmax, -vmax/4., nz_1+1, endpoint=True)[1:]
nz_2 = 30
vz_2 = np.linspace(-vmax/4., vmax/4., nz_2+1, endpoint=True)[1:]
nz_3 = 15
vz_3 = np.linspace(vmax/4., vmax, nz_3+1, endpoint=True)[1:]
nz = nz_1 + nz_2 + nz_3
vz = np.r_[vz_1, vz_2, vz_3]
nz_array = np.arange(nz)

x_out = vx * texp
y_out = vy * texp
z_out = vz * texp

rho  = np.zeros((nx,ny,nz))
temp = np.zeros((nx,ny,nz))
comp = np.zeros((nx,ny,nz,len(Z)))
erad = np.zeros((nx,ny,nz))
vxx  = np.zeros((nx,ny,nz))
vxy  = np.zeros((nx,ny,nz))
vxz  = np.zeros((nx,ny,nz))

for i in range(nx):
	for j in range(ny):
		for k in range(nz):
			vr = (vx[i]**2 + vy[j]**2 + vz[k]**2)**0.5

			if (vr < vmax): 
				rho[i,j,k] = rho0
				temp[i,j,k] = T0
			else:
				rho[i,j,k] = rho0*1e-20
				temp[i,j,k] = T0*1e-4

			vxx[i,j,k] = vx[i];
			vxy[i,j,k] = vy[j];
			vxz[i,j,k] = vz[k];

			# get composition
			m_enc = 4.0*pi/3.0*(vr*texp)**3.0*rho0/m_sun	
			if   (m_enc < 0.50): ni_frac = 1.0
			elif (m_enc < 0.75): ni_frac = (0.75 - m_enc)/0.25
			else: ni_frac = 0
			comp[i][j][k][0] = 1 - ni_frac
			comp[i][j][k][1] = 0.0
			comp[i][j][k][2] = 0.0
			comp[i][j][k][3] = ni_frac;

fout = h5py.File('lucy_3D-testing_nonuniform_grid.h5','w')
fout.create_dataset('time',data=[texp],dtype='d')
fout.create_dataset('Z',data=Z,dtype='i')
fout.create_dataset('A',data=A,dtype='i')
fout.create_dataset('rho',data=rho,dtype='d')
fout.create_dataset('temp',data=temp,dtype='d')
fout.create_dataset('vx',data=vxx,dtype='d')
fout.create_dataset('vy',data=vxy,dtype='d')
fout.create_dataset('vz',data=vxz,dtype='d')
fout.create_dataset('erad',data=erad,dtype='d')
fout.create_dataset('comp',data=comp,dtype='d')
fout.create_dataset('x_out',data=x_out,dtype='d')
fout.create_dataset('y_out',data=y_out,dtype='d')
fout.create_dataset('z_out',data=z_out,dtype='d')
fout.create_dataset('rmin',data=[-1.0*rmax,-1.0*rmax,-1.0*rmax],dtype='d')