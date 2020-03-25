import h5py
import pylab as py

#############################
## make a simple atom
#############################

fname = 'test_atom.hdf5'

# levels
level_g     = [1,    1,      1,        1,   1]
level_E     = [0,  1.0,   10.0,     11.0, 0.0]
level_i     = [0,    0,      0 ,       0,   1]

#ionization energies
ion_chi     = [20.0, 999999]
# index of ground states of all ionization stages
ion_ground  = [0, 4]

# lines
lu          = [  1, 3 ,2]
ll          = [  0 ,2, 0]
Aij         = [1e9,1e9,1]


f = h5py.File(fname,'w')

base = "1/"
grp = f.create_group(base)

f[base].attrs["n_ions"]    = len(ion_chi)
f[base].attrs["n_levels" ] = int(len(level_g))
f[base].attrs["n_lines" ]  = len(Aij)

f.create_dataset(base + "ion_chi",data = ion_chi)
f.create_dataset(base + "ion_ground",data = ion_ground)
f.create_dataset(base + "level_g",data = level_g)
f.create_dataset(base + "level_E",data = level_E)
f.create_dataset(base + "level_i",data = level_i)
f.create_dataset(base + "line_l",data = ll)
f.create_dataset(base + "line_u",data = lu)
f.create_dataset(base + "line_A",data = Aij)
f.close()

print ("lines:")
for i in range(len(lu)):
    h_ev = 4.135667696e-15
    dE = level_E[lu[i]] - level_E[ll[i]]
    nu = dE/h_ev
    print lu[i],ll[i],dE,nu,Aij[i]
