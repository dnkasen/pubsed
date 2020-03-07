import numpy as np
import scipy as sp
import h5py

##### Constants
m_sun  = 1.99e33 # [g]
c = 2.99e10      # [cm/s]

##### Conversion factors
dtos = 24.*3600.    # multiply by this to convert from days to seconds
htos = 3600.        # multiply by this to convert from hours to seconds

##### Basic simulation parameters
nx = 100 # number of zones in the x-direction
nz = 200 # number of zones in the z-direction
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

dvx = ax/(1.0*nx)
dvz = (2.0*az)/(1.0*nz)

nx_array = np.arange(nx)
nz_array = np.arange(nz)

vx = np.array([(1 + i) * dvx for i in np.arange(nx)])
vz = np.array([(1 + k - nz/2) * dvz for k in np.arange(nz)])
vxmid = np.array([(0.5 + i) * dvx for i in np.arange(nx)])
vzmid = np.array([(0.5 + k - nz/2) * dvz for k in np.arange(nz)])

rho  = np.zeros((nx,nz))
temp = np.zeros((nx,nz))
comp = np.zeros((nx,nz,len(Zlist)))
erad = np.zeros((nx,nz))
velx = np.zeros((nx,nz))
velz = np.zeros((nx,nz))



### Fill the arrays

M_calc = 0
KE_calc = 0
V0_calc = 0

for i in nx_array:
    for k in nz_array:

        velx[i,k] = vx[i];
        velz[i,k] = vz[k];

        vrmidsquared = vxmid[i]**2. + vzmid[k]**2.
        ssquared = (vxmid[i]/ax)**2. + (vzmid[k]/az)**2.

        if ssquared < 1:

            rho[i,k] = rho0*((ssquared**(1./2.))**alpha)
            temp[i,k] = T0

            if i == 0:
                M_calc += rho[i,k]*((np.pi*((vx[i]*texp)**2.))*(dvz*texp))
                KE_calc += (1./2.)*rho[i,k]*((np.pi*((vx[i]*texp)**2.))*(dvz*texp))*vrmidsquared
                V0_calc += (np.pi*((vx[i]*texp)**2.))*(dvz*texp)
            else:
                M_calc += rho[i,k]*((np.pi*((vx[i]*texp)**2.-(vx[i-1]*texp)**2.))*(dvz*texp))
                KE_calc += (1./2.)*rho[i,k]*((np.pi*((vx[i]*texp)**2.-(vx[i-1]*texp)**2.))*(dvz*texp))*vrmidsquared
                V0_calc += (np.pi*((vx[i]*texp)**2. - (vx[i-1]*texp)**2.))*(dvz*texp)

        for n in range(nelems):
            comp[i,k,n] = comp_uniform[n]



### Write the model file

fout = h5py.File(name + '_2D.h5','w')
fout.create_dataset('time',data=[texp],dtype='d')
fout.create_dataset('Z',data=Zlist,dtype='i')
fout.create_dataset('A',data=Alist,dtype='i')
fout.create_dataset('rho',data=rho,dtype='d')
fout.create_dataset('temp',data=temp,dtype='d')
fout.create_dataset('vx',data=velx,dtype='d')
fout.create_dataset('vz',data=velz,dtype='d')
fout.create_dataset('erad',data=erad,dtype='d')
fout.create_dataset('comp',data=comp,dtype='d')
fout.create_dataset('dr',data=[dvx*texp,dvz*texp],dtype='d')