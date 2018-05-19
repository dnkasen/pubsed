import numpy as np
import h5py

####################################
m_sun  = 1.99e33
pi     = 3.14159
nx     = 100
mass   = 1.0*1.99e33
vmax   = 1.0e9
texp   = 20.0*(3600.0*24.0)
rmax  = vmax*texp
Z  = [1]
A  = [1]
T0     = 1.0e4*(20)
rho0    = 1e-20
###################################

##################################
# Make sedona 1D spherical model
#################################

fout = open("vacuum_1D.mod","w")

dv    = vmax/(1.0*nx)
rmin = 0.0
fout.write("1D_sphere standard\n")
fout.write(str(nx) + "\t" + str(rmin) + "\t" + str(texp) + " 1 \n")
fout.write("1.1\n")

for i  in range(nx):
    v = (i+1.0)*dv
    r = v*texp
    line = "%10.4e %10.4e %10.4e %10.4e 1.0\n" % (r,0,rho0,T0)
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
		vr = vx[i]**2 + vz[j]**2
		if (vr < vmax*vmax): 
			rho[i,j] = rho0
			temp[i,j] = T0
		else:
			rho[i,j] = rho0*1e-20
			temp[i,j] = T0*1e-4

		vxx[i,j] = vx[i];
		vxz[i,j] = vz[j];
		comp[i,j,0] = 1
	

fout = h5py.File('vacuum_2d.h5','w')
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
			if (vr < vmax*vmax): 
				rho[i,j,k] = rho0
				temp[i,j,k] = T0
			else:
				rho[i,j,k] = rho0*1e-20
				temp[i,j,k] = T0*1e-4

			vxx[i,j,k] = vx[i];
			vxy[i,j,k] = vy[j];
			vxz[i,j,k] = vz[j];
			comp[i,j,k,0] = 1

fout = h5py.File('vacuum_3d.h5','w')
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

