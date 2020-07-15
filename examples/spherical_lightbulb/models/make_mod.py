import numpy as np
import h5py

####################################
nx     = 100
texp   = 1
rmax  = 3e15
rph   = 6e14
Z  = [1]
A  = [1]
T0     = 1.0e4
rho0   = 1.0/(rmax-rph)
###################################

##################################
# Make sedona 1D spherical model
#################################

fout = open("vacuum_1D.mod","w")

dr    = rmax/(1.0*nx)
fout.write("1D_sphere standard\n")
fout.write(str(nx) + "\t" + str(0) + "\t" + str(texp) + " 1 \n")
fout.write("1.1\n")

for i  in range(nx):
    v = 0.0
    r = dr*(i+1)
    line = "%10.4e %10.4e %10.4e %10.4e 1.0\n" % (r,0,rho0,T0)
    fout.write(line)

fout.close()


##################################
# Make sedona 2D hdf5 model
#################################

dr   = rmax/(1.0*nx)
nz   = nx*2
vx   = np.zeros((nx,nz))
vz   = np.zeros((nx,nz))
rho  = np.zeros((nx,nz))
temp = np.zeros((nx,nz))
comp = np.zeros((nx,nz,len(Z)))
erad = np.zeros((nx,nz))
vxz  = np.zeros((nx,nz))
vxx  = np.zeros((nx,nz))

for i in range(nx):
    for j in range(nz):

        x = dr*(i+0.5)
        z = dr*(j+0.5) - rmax
        rsq = x**2 + z**2
        if (rsq < rmax*rmax):
            rho[i,j] = rho0
            temp[i,j] = T0
        else:
            rho[i,j] = rho0*1e-20
            temp[i,j] = T0*1e-4

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
fout.create_dataset('dr',data=[dr,dr],dtype='d')



#########################################
# Make sedona 3D Cartesian hdf5 model
#########################################
dr   = rmax/(1.0*nx)*2

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

            x = (i+0.5) - rmax
            y = (j+0.5) - rmax
            z = (k+0.5) - rmax
            rsq = x**2 + y**2 + z**2
            if (rsq < rmax*rmax):
                rho[i,j,k] = rho0
                temp[i,j,k] = T0
            else:
                rho[i,j,k] = rho0*1e-20
                temp[i,j,k] = T0*1e-4


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
fout.create_dataset('dr',data=[dr,dr,dr],dtype='d')
fout.create_dataset('rmin',data=[-1.0*rmax,-1.0*rmax,-1.0*rmax],dtype='d')
