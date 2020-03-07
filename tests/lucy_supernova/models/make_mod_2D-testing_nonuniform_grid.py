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
# Make sedona 2D hdf5 model
#################################
nx_1 = 25
vx_1 = np.linspace(0, vmax/4., nx_1+1, endpoint=True)[1:]
nx_2 = 25
vx_2 = np.linspace(vmax/4., vmax, nx_2+1, endpoint=True)[1:]
nx = nx_1 + nx_2
vx = np.r_[vx_1, vx_2]
nx_array = np.arange(nx)

nz_1 = 25
vz_1 = np.linspace(-vmax, -vmax/4., nz_1+1, endpoint=True)[1:]
nz_2 = 50
vz_2 = np.linspace(-vmax/4., vmax/4., nz_2+1, endpoint=True)[1:]
nz_3 = 25
vz_3 = np.linspace(vmax/4., vmax, nz_3+1, endpoint=True)[1:]
nz = nz_1 + nz_2 + nz_3
vz = np.r_[vz_1, vz_2, vz_3]
nz_array = np.arange(nz)

x_out = vx * texp
z_out = vz * texp

rho  = np.zeros((nx,nz))
temp = np.zeros((nx,nz))
comp = np.zeros((nx,nz,len(Z)))
erad = np.zeros((nx,nz))
vxz  = np.zeros((nx,nz))
vxx  = np.zeros((nx,nz))

for i in range(nx):
    for j in range(nz):
        vr = (vx[i]**2 + vz[j]**2)**0.5
        if (vr < vmax):
            rho[i,j] = rho0
            temp[i,j] = T0
        else:
            rho[i,j] = rho0*1e-20
            temp[i,j] = T0*1e-4

        vxx[i,j] = vx[i];
        vxz[i,j] = vz[j];

        # get composition
        m_enc = 4.0*pi/3.0*(vr*texp)**3.0*rho0/m_sun
        if   (m_enc < 0.50): ni_frac = 1.0
        elif (m_enc < 0.75): ni_frac = (0.75 - m_enc)/0.25
        else: ni_frac = 0
        comp[i][j][0] = 1 - ni_frac
        comp[i][j][1] = 0.0
        comp[i][j][2] = 0.0
        comp[i][j][3] = ni_frac;


fout = h5py.File('lucy_2D-testing_nonuniform_grid.h5','w')
fout.create_dataset('Z',data=Z,dtype='i')
fout.create_dataset('A',data=A,dtype='i')
fout.create_dataset('rho',data=rho,dtype='d')
fout.create_dataset('temp',data=temp,dtype='d')
fout.create_dataset('vx',data=vxx,dtype='d')
fout.create_dataset('vz',data=vxz,dtype='d')
fout.create_dataset('erad',data=temp,dtype='d')
fout.create_dataset('comp',data=comp,shape=comp.shape,dtype='f')
fout.create_dataset('time',data=[texp],dtype='d')
fout.create_dataset('x_out',data=x_out,dtype='d')
fout.create_dataset('z_out',data=z_out,dtype='d')
fout.create_dataset('rmin',data=[0.0,-1.0*rmax],dtype='d')