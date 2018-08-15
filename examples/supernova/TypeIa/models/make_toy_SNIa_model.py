import numpy as np
import h5py

#######################################################
nx     = 100         # number of radial zones
n_3d   =  60         # number zones per side for 3D model
mass   = 1.0         # mass in solar masses
KE     = 1.0e51      # kinetic energy in ergs
texp   = 2.0         # time since explosion in days
vmax   = 3.0e9       # outer velocity (cm/s)
Mni    = 0.6         # mass of 56Ni
Mime   = 0.3         # mass of IME
T0     = 1e4         # initial temperatures
name   = "toy_SNIa"  # base name of model
#######################################################

# elements included
Z = [ 6, 8,14,16,20,26,27,28]
A = [12,16,28,32,40,56,56,56]

comp_ni  = [  0.,  0.,   0.,   0.,   0.,    0., 0,  1.]
comp_ime = [  0.,  0., 0.7, 0.29,  0.01,    0., 0., 0.]
comp_co  = [0.5, 0.5,   0.,    0., 1e-4,  1e-3, 0., 0.]

# pre calculations
pi     = 3.14159
m_sun  = 1.99e33
rmax  = vmax*texp
v_e   = (KE/6/(mass*m_sun))**(0.5)
rho0  = mass*m_sun/8/pi/(v_e*texp)**3.0
nelems = len(Z)
dv    = vmax/(1.0*nx)

v1d     = np.arange(dv,vmax*1.0001,dv)
t1d     = np.zeros(nx)
rho1d   = np.zeros(nx)
comp_1d = np.zeros((nx,nelems))

xx = v1d/v_e
mfrac = 1 - 0.5*np.exp(-xx)*(xx**2 + 2*xx + 2)
mmid = 0.5*(Mni + Mime + Mni)/mass
mdel = Mime/mass/1.85
ximem = np.exp(-1.0*(mfrac - mmid)**6.0/mdel**6.0)


##################################
# Make sedona 1D spherical model
#################################

fout = open(name + "_1D.mod","w")

rmin = 0.0
fout.write("1D_sphere standard\n")
fout.write(str(nx) + "\t" + str(rmin) + "\t" + str(texp) + " " + str(len(Z)) + " \n")

# write out elements
for xZ, xA in zip(Z,A):
	fout.write(str(xZ) + "." + str(xA) + " "),

# write out layer properties
for i  in range(nx):
    v  = (i+1.0)*dv
    vm = (i+0.5)*dv
    r = v*texp
    rho1d[i] = rho0*np.exp(-vm/v_e)

    if (vm < 1.0e9):
    	t1d[i] = T0*(1 + (1e9 - vm)/1e9)
    else:
    	t1d[i] = T0*(vm/1e9)**(-2.0)
    
    line = "%10.4e %10.4e %10.4e %10.4e " % (r,v,rho1d[i],t1d[i])
    
    xime = ximem[i]
    xni = 0
    xco = 0
    if (mfrac[i] < mmid): 
    	xni = 1 - xime
    else:  
		xco = 1 - xime

    for j in range(nelems):
    	comp_1d[i,j] = xni*comp_ni[j] + xime*comp_ime[j] + xco*comp_co[j]
    	line = line + str(comp_1d[i,j]) + " "
    fout.write("\n")
    fout.write(line)

fout.close()

##################################
# Make sedona 2D hdf5 model
#################################

dv    = vmax/(1.0*nx)
vx   = np.arange(dv,vmax+0.1,dv)
vz   = np.arange(-1.0*vmax + dv,vmax+0.1,dv)
nz   = nx*2
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
			ind = np.searchsorted(v1d,vr)
			comp[i,j,:] = comp_1d[ind,:]
			temp[i,j]   = t1d[ind]
			rho[i,j]    = rho1d[ind]


		else:
			rho[i,j]  = rho0*1e-10
			temp[i,j] = T0*1e-4
			comp[i,j,:] = comp_co


		vxx[i,j] = vx[i];
		vxz[i,j] = vz[j];
	

fout = h5py.File(name + '_2D.h5','w')
fout.create_dataset('time',data=[texp],dtype='d')
fout.create_dataset('Z',data=Z,dtype='i')
fout.create_dataset('A',data=A,dtype='i')
fout.create_dataset('rho',data=rho,dtype='d')
fout.create_dataset('temp',data=temp,dtype='d')
fout.create_dataset('vx',data=vxx,dtype='d')
fout.create_dataset('vz',data=vxz,dtype='d')
fout.create_dataset('erad',data=erad,dtype='d')
fout.create_dataset('comp',data=comp,dtype='d')
fout.create_dataset('dr',data=[dv*texp,dv*texp],dtype='d')


#########################################
# Make sedona 3D Cartesian hdf5 model
#########################################
dv    = vmax/(1.0*nx)*2.0
vx    = np.arange(-1.0*vmax + dv,vmax+0.1,dv)
vy    = np.arange(-1.0*vmax + dv,vmax+0.1,dv)
vz    = np.arange(-1.0*vmax + dv,vmax+0.1,dv)

nx = n_3d
ny = nx
nz = nx
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
			vr = vx[i]**2 + vy[j]**2 + vz[k]**2

			if (vr < vmax): 
				ind = np.searchsorted(v1d,vr)
				comp[i,j,:] = comp_1d[ind,:]
				temp[i,j]   = t1d[ind]
				rho[i,j]    = rho1d[ind]

			else:
				rho[i,j]  = rho0*1e-10
				temp[i,j] = T0*1e-4
				comp[i,j,:] = comp_co

				vxx[i,j,k] = vx[i];
				vxy[i,j,k] = vy[j];
				vxz[i,j,k] = vz[j];
				comp[i,j,k,0] = 1

fout = h5py.File(name + '_3D.h5','w')
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
fout.create_dataset('dr',data=[dv*texp,dv*texp,dv*texp],dtype='d')
fout.create_dataset('rmin',data=[-1.0*rmax,-1.0*rmax,-1.0*rmax],dtype='d')

