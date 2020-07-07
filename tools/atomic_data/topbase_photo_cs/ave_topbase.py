import numpy as np
import matplotlib.pyplot as plt
import h5py

rydberg = 13.605693122994

fname = "topbase_HI_photocs.dat"
fin = open(fname,"r")

# read off header
fin.readline()
fin.readline()
fin.readline()

fout = h5py.File("topbase_PI.hdf5","w")

base = "1/0"
fout.create_group(base)

# read data until we hit last line with a "=" in it
data = (fin.readline()).split()

nval = 1
ncnt = 0
cnt = 0
norm = 0

# empty arrays to hold cross-section
E_cs = np.zeros(shape=(0),dtype='d')
s_cs = np.zeros(shape=(0),dtype='d')

while (not "=" in data[0]):

    # threshold energy to ionize level
    E_thr = -1*(float(data[5]))*rydberg
    # level energy relative to ground
    Elev = rydberg - E_thr
    # level P value
    lval = int((data[3])[1])
    ilev = int((data[4]))

    gval = 2*(2*lval + 1)
    # read cross-section
    i = 0
    while (True):
        data = (fin.readline()).split()

        # if not two element is line, have reached end
        if (len(data) != 2):
            break

        this_E = float(data[0])*rydberg
        this_s = gval*float(data[1])*1e-18

        # add into average
        if (ncnt == 0):
            E_cs = np.append(E_cs,this_E)
            s_cs = np.append(s_cs,this_s)
        else:
            s_cs[i] += this_s

        i += 1

    ncnt += 1
    norm += gval
    # done wiht this n keep it going
    if (ncnt == nval):
        # average
        s_cs /= (1.0*norm)

        print nval,norm

        name = base + "/cs_" + str(cnt)
        fout.create_group(name)
        fout.create_dataset(name + "/E_ev",data=E_cs)
        fout.create_dataset(name + "/sigma",data=s_cs)
        fout.create_dataset(name + "/n_pts",data=len(E_cs))
        nval += 1
        ncnt = 0
        cnt += 1
        norm = 0

        # empty arrays to hold cross-section
        E_cs = np.zeros(shape=(0),dtype='d')
        s_cs = np.zeros(shape=(0),dtype='d')





fout.create_dataset(base + "/n_cs",data=cnt)

exit(0)
