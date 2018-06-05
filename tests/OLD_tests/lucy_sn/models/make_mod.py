import numpy as np
import h5py

####################################
m_sun  = 1.99e33
pi     = 3.14159
nx     = 100
mass   = 1.4*1.99e33
vmax   = 1.0e9
texp   = 1.0*(3600.0*24.0)
elems  = ('14.28','26.56','27.56','28.56')
T0     = 0.0
###################################


rmax    = vmax*texp
rho0    = mass/(4.0*pi/3.0*rmax**3)
rho_min = rho0*1e-10
n_elems = len(elems)


##################################
# Make 1D hdf5 model
##################################

dv      = vmax/(1.0*nx)

fout = open("lucy_1D.mod","w")
fout.write("1D_sphere SNR\n")
fout.write(str(nx) + "\t" + str(0.0) + "\t" + str(texp) + "\t")
fout.write(str(n_elems) + "\n")
for el in elems: fout.write((el) + " ")
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
exit(0)


##################################
# Make 3D hdf5 model
##################################

dv      = vmax*2.0/(1.0*nx)


# open hdf5 file
f = h5py.File('lucy_3D.hdf5','w')

# write attributes
f.attrs.create("texp",texp)
f.attrs.create("dvel",dv)
f.attrs.create("n_x",nx,dtype='i')
f.attrs.create("n_elems",n_elems,dtype='i')

# write element list
f.create_dataset('elems',data=elems,dtype='i')

# set up model density, temp, composition
comp = np.zeros((nx,nx,nx,n_elems))
dens = np.zeros((nx,nx,nx))
temp = np.zeros((nx,nx,nx))
count = 0.0
for i in range(nx):
    for j in range(nx):
        for k in range(nx):

            vx = i*dv - vmax
            vy = j*dv - vmax
            vz = k*dv - vmax
            vr = np.sqrt(vx**2 + vy**2 + vz**2)

            if (vr < vmax): dens[i,j,k] = rho0
            else:           dens[i,j,k] = rho_min

            m_enc = 4.0*pi/3.0*(vr*texp)**3.0*rho0/m_sun
            if   (m_enc < 0.50): ni_frac = 1.0
            elif (m_enc < 0.75): ni_frac = (0.75 - m_enc)/0.25
            else: ni_frac = 0

            comp[i,j,k,0] = 1 - ni_frac;
            comp[i,j,k,1] = 0
            comp[i,j,k,2] = 0
            comp[i,j,k,0] = ni_frac;
            
            temp[i,j,k] = T0


# write out data

dset = f.create_dataset("dens",data = dens)
dset = f.create_dataset("temp",data = temp)
for i in range(n_elems):
    name = "comp" + str(elems[i])
    dset = f.create_dataset(name,data=comp[:,:,:,i])
f.close()
