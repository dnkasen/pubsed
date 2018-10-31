import numpy as np
import h5py

####################################
m_sun  = 1.99e33
pi     = 3.14159
n_1d   = 100
n_2d   = 70
n_3d   = 50
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
# Make 1D ascii model
##################################
nx = n_1d
dv      = vmax/(1.0*nx)
fout = open("lucy_1D.mod","w")
fout.write("1D_sphere SNR\n")
fout.write(str(nx) + "\t" + str(0.0) + "\t" + str(texp) + "\t")
fout.write(str(n_elems) + "\n")
for this_Z,this_A in zip(Z,A): fout.write(str(this_Z) + "." + str(this_A) + " ")
fout.write("\n")

comp = np.zeros(4)
for i in range(nx):
    v = dv*(i+1.0)
    line = "%10.4e %10.4e %10.4e " % (v,rho0,T0)
    fout.write(line)

    # get composition
    m_enc = 4.0*pi/3.0*(v*texp)**3.0*rho0/m_sun
    if   (m_enc < 0.50): ni_frac = 1.0
    elif (m_enc < 0.75): ni_frac = (0.75 - m_enc)/0.25
    else: ni_frac = 0
    comp[0] = 1 - ni_frac;
    comp[1] = 0
    comp[2] = 0
    comp[3] = ni_frac;

    # write composition
    for c in comp: fout.write("%10.4e " % c)
    fout.write("\n")

fout.close()

##################################
# Make sedona 2D hdf5 model
#################################
nx = n_2d
nz = n_2d*2
dv    = vmax/(1.0*nx)
vx  = np.arange(dv,vmax+0.1,dv)
vz  = np.arange(-1.0*vmax + dv,vmax+0.1,dv)
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


fout = h5py.File('lucy_2D.h5','w')
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
