import numpy as np
import h5py

####################################
m_sun  = 1.99e33
pi     = 3.14159
nx     = 400
nx_2d  = 100
mass   = 1.0*1.99e33
vmax   = 1.0e9
texp   = 1.0*(3600.0*24.0)
elems  = ('26.56','27.56','28.56')
Z  = [26,27,28]
A  = [56,56,56]
T0     = 1.0e4*(20)
###################################


rmax    = vmax*texp
rho0    = mass/(4.0*pi/3.0*rmax**3)
n_elems = len(elems)


##################################
# Make sedona6 1D ascii model
##################################

dv    = vmax/(1.0*nx)
v     = np.arange(dv,vmax+0.1,dv)
rho   = 0*v + rho0
temp  = 0*v + T0
comp0 = 0*v 
comp1 = 0*v
comp2 = 0*v + 1.0

fout = open("nickelball_1D.mod","w")
fout.write("1D_sphere SNR\n")
fout.write(str(nx) + "\t" + str(0.0) + "\t" + str(texp) + "\t")
fout.write(str(n_elems) + "\n")
for el in elems: fout.write((el) + " ")
fout.write("\n")

for i in range(nx):
    line = "%10.4e %10.4e %10.4e " % (v[i],rho[i],temp[i])
    fout.write(line)
    line = "%10.4e %10.4e %10.4e\n" % (comp0[i],comp1[i],comp2[i])
    fout.write(line)

fout.close()


##################################
# Make  1D sedonabox hdf5 model
##################################


# openfile
fout = h5py.File('nickelball_1d.h5','w')
fout.create_dataset('Version',data=[0],dtype='i')
fout.create_dataset('Z',data=Z,dtype='i')
fout.create_dataset('A',data=A,dtype='i')
fout.create_dataset('vx',data=v,dtype='d')
fout.create_dataset('rho',data=rho,dtype='d')
fout.create_dataset('time',data=texp,dtype='d')
fout.create_dataset('comp0',data=comp0,dtype='d')
fout.create_dataset('comp1',data=comp1,dtype='d')
fout.create_dataset('comp2',data=comp2,dtype='d')



##################################
# Make sedona 2D hdf5 model
#################################
nx = nx_2d
dv    = vmax/(1.0*nx)
vx  = np.arange(dv,vmax+0.1,dv)
vz  = np.arange(-1.0*vmax + dv,vmax+0.1,dv)
nx = nx
nz = nx*2
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
		if (m_enc < 0.50): ni_frac = 1.0
		elif (m_enc < 0.75): ni_frac = (0.75 - m_enc)/0.25
		else: ni_frac = 0
		comp[i][j][0] = 0.0
		comp[i][j][1] = 0.0
		comp[i][j][2] = 1.0



fout = h5py.File('nickelball_2d.h5','w')
fout.create_dataset('Z',data=Z,dtype='i')
fout.create_dataset('A',data=A,dtype='i')
fout.create_dataset('rho',data=rho,dtype='d')
fout.create_dataset('temp',data=temp,dtype='d')
fout.create_dataset('vx',data=vxx,dtype='d')
fout.create_dataset('vz',data=vxz,dtype='d')
fout.create_dataset('erad',data=temp,dtype='d')
fout.create_dataset('comp',data=comp,shape=comp.shape,dtype='f')
fout.create_dataset('dr',data=[dv*texp,dv*texp],dtype='d')
fout.create_dataset('time',data=[texp],dtype='d')


