import numpy as np
import scipy as sp
import h5py
import matplotlib.pyplot as plt

##### Constants
m_sun  = 1.99e33 # [g]
c = 2.99e10      # [cm/s]

##### Conversion factors
dtos = 24.*3600.    # multiply by this to convert from days to seconds
htos = 3600.        # multiply by this to convert from hours to seconds

##### Basic simulation parameters
# nx = 100 # number of zones in the x-direction
# nz = 200 # number of zones in the z-direction
T0 = 1.0e4 # initial temperatures [K]
texp = 0.1*dtos # [s]
name = "ellipsoidal_outflow"  # base name of model

### Elements included
Zlist = [26,27,28]
Alist = [56,56,56]
nelems = len(Zlist)

### Composition
comp_uniform = [0.,0.,1.]

### Simulation parameters
M = (1.e-1)*m_sun                   # mass in [g]
vchar = (1.e-1)*c           # characteristic velocity of the ejecta
axial_ratio = 4.            # axial ratio, i.e. the ratio of ax/az, assuming the ax > az
KE = (1./2.)*M*(vchar**2.)  # kinetic energy [ergs]

### Density profile input parameters
alpha = -1. # exponent in the power law

### Ellipsoid quantities. The ellipsoid is an oblate spheroid
# ax, az are the semi-major axes of the ellipsoid
az = ((6.*KE/M)/(1 + 2.*(axial_ratio**2.)))**(1./2.) * ( (1./(3.+alpha)) / (1./(5.+alpha)) )**(1./2.)
ax = az*axial_ratio
vmax = np.amax(np.array([ax, az]))
xmax0 = ax*texp
zmax0 = az*texp
rmax0 = vmax*texp

### Global quantities
V0 = (4./3.)*np.pi*(xmax0**2.)*zmax0    # the volume [cm^3] at t=texp
rho0 = M / ( 4.*np.pi*(((ax**2.)*az)*(texp**3.)) * (1./(3.+alpha)) ) # density scale [g/cm^3] at t=texp



##### Write the outputs

### Set the array sizes

# dvx = ax/(1.0*nx)
# dvz = (2.0*az)/(1.0*nz)

# nx_array = np.arange(nx)
# nz_array = np.arange(nz)

nx_1 = 25
vx_1 = np.linspace(0, ax/4., nx_1+1, endpoint=True)[1:]
nx_2 = 25
vx_2 = np.linspace(ax/4., ax, nx_2+1, endpoint=True)[1:]
nx = nx_1 + nx_2
vx = np.r_[vx_1, vx_2]
nx_array = np.arange(nx)

nz_1 = 25
vz_1 = np.linspace(-az, -az/4., nz_1+1, endpoint=True)[1:]
nz_2 = 50
vz_2 = np.linspace(-az/4., az/4., nz_2+1, endpoint=True)[1:]
nz_3 = 25
vz_3 = np.linspace(az/4., az, nz_3+1, endpoint=True)[1:]
nz = nz_1 + nz_2 + nz_3
vz = np.r_[vz_1, vz_2, vz_3]
nz_array = np.arange(nz)

vxmid = np.r_[[0], vx]
vxmid = ((vxmid + np.roll(vxmid, -1))/2.)[:-1]

vzmid = np.r_[[-az], vz]
vzmid = ((vzmid + np.roll(vzmid, -1))/2.)[:-1]

x_out = vx * texp
z_out = vz * texp

x_edges = np.r_[[0], x_out]
z_edges = np.r_[[-zmax0], z_out]

rho  = np.zeros((nx,nz))
temp = np.zeros((nx,nz))
comp = np.zeros((nx,nz,len(Zlist)))
erad = np.zeros((nx,nz))
velx = np.zeros((nx,nz))
velz = np.zeros((nx,nz))
velxmid = np.zeros((nx,nz))
velzmid = np.zeros((nx,nz))



### Fill the arrays

M_calc = 0
KE_calc = 0
V0_calc = 0

for i in nx_array:
    for k in nz_array:

        velx[i,k] = vx[i];
        velz[i,k] = vz[k];

        velxmid[i,k] = vxmid[i];
        velzmid[i,k] = vzmid[k];

        vrmidsquared = vxmid[i]**2. + vzmid[k]**2.
        ssquared = (vxmid[i]/ax)**2. + (vzmid[k]/az)**2.

        if ssquared < 1:

            rho[i,k] = rho0*((ssquared**(1./2.))**alpha)
            temp[i,k] = T0

            M_calc += rho[i,k]*((np.pi*(x_edges[i+1]**2. - x_edges[i]**2.))*(z_edges[k+1]-z_edges[k]))
            KE_calc += (1./2.)*rho[i,k]*((np.pi*(x_edges[i+1]**2. - x_edges[i]**2.))*(z_edges[k+1]-z_edges[k]))*vrmidsquared
            V0_calc += ((np.pi*(x_edges[i+1]**2. - x_edges[i]**2.))*(z_edges[k+1]-z_edges[k]))

        for n in range(nelems):
            comp[i,k,n] = comp_uniform[n]

# print(M, M_calc)
# print(KE, KE_calc)
# print(V0, V0_calc)



### Write the model file

fout = h5py.File(name + '_2D-testing_nonuniform_grid.h5','w')
fout.create_dataset('time',data=[texp],dtype='d')
fout.create_dataset('Z',data=Zlist,dtype='i')
fout.create_dataset('A',data=Alist,dtype='i')
fout.create_dataset('rho',data=rho,dtype='d')
fout.create_dataset('temp',data=temp,dtype='d')
fout.create_dataset('vx',data=velx,dtype='d')
fout.create_dataset('vz',data=velz,dtype='d')
fout.create_dataset('erad',data=erad,dtype='d')
fout.create_dataset('comp',data=comp,dtype='d')
fout.create_dataset('x_out',data=x_out,dtype='d')
fout.create_dataset('z_out',data=z_out,dtype='d')
fout.create_dataset('rmin',data=[0.0,-1.0*zmax0],dtype='d')



### Plot the density vs velocity

# x-direction
plt.plot(velxmid[:,int(nz/2)], rho[:,int(nz/2)], 'k')
plt.plot(velxmid[:,int(nz/2)], rho[:,int(nz/2)], 'ko', markersize=2)
for vc in velx[:,int(nz/2)]:
    plt.axvline(x=vc, color='k', linestyle='--', linewidth=0.5)
plt.xlabel('$v_x$ [cm/s] (roughly along $z=0$ axes)')
plt.ylabel('$\\hat{\\rho}(v_x,t_{exp})$ [g/cm$^3$]')
plt.savefig('rho_vs_vx-testing_nonuniform_grid.pdf')
plt.close()

# z-direction
plt.plot(velzmid[0,:], rho[0,:], 'k')
plt.plot(velzmid[0,:], rho[0,:], 'ko', markersize=2)
for vc in velz[0,:]:
    plt.axvline(x=vc, color='k', linestyle='--', linewidth=0.5)
plt.xlabel('$v_z$ [cm/s] (roughly along the $x = y = 0$)')
plt.ylabel('$\\hat{\\rho}(v_z,t_{exp})$ [g/cm$^3$]')
plt.savefig('rho_vs_vz-testing_nonuniform_grid.pdf')
plt.close()