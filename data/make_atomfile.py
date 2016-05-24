import h5py
import pylab as py


f = h5py.File('2level_atomdata.hdf5','w')

base = "1/"
grp = f.create_group(base)

g = [2,   6, 1]
E = [0,10.2, 0]
i = [0,   0, 1]

ll = [0]
lu = [1]
Aij = [1e9] 

chi = [13.6,99999]
gnd = [0,2]
f[base].attrs["n_ions"]    = len(chi)
f[base].attrs["n_levels" ] = int(len(g))
f[base].attrs["n_lines" ]  = len(Aij)


f.create_dataset(base + "ion_chi",data = chi)
f.create_dataset(base + "ion_ground",data = gnd)
f.create_dataset(base + "level_g",data = g)
f.create_dataset(base + "level_E",data = E)
f.create_dataset(base + "level_i",data = i)
f.create_dataset(base + "line_l",data = ll)
f.create_dataset(base + "line_u",data = lu)
f.create_dataset(base + "line_A",data = Aij)

f.close()

g = h5py.File('atomdata.hdf5','r')
print g["1/"].attrs["n_levels"]
