import h5py
import numpy as np

fout = h5py.File("gaunt_ff_data.hdf5","w")

for Z in range(1,37):

    name = "gauntff_merged_Z0" + str(Z) + ".dat"
    if (Z >= 10):
        name = "gauntff_merged_Z" + str(Z) + ".dat"

    print name
    fin = open(name,"r")
    line = fin.readline()
    line = fin.readline()
    data = line.split()
    n_g = int(data[0])
    n_u   = int(data[1])
    data = (fin.readline()).split()
    log_g_start = float(data[0])
    data = (fin.readline()).split()
    log_u_start = float(data[0])
    data = (fin.readline()).split()
    d_gu = float(data[0])
    fin.readline()
    fin.readline()
    fin.readline()

    gaunt = np.zeros((n_u,n_g))
    for i_u in range(0,n_u):
        data = (fin.readline()).split()
        gaunt[i_u,:] = data

    log_g = log_g_start + d_gu*np.arange(0,n_g)
    log_u = log_u_start + d_gu*np.arange(0,n_u)

    fout.create_group(str(Z))
    fout.create_dataset(str(Z) + "/log_gam",data=log_g)
    fout.create_dataset(str(Z) + "/log_u",data=log_u)
    fout.create_dataset(str(Z) + "/gaunt",data=gaunt)



fout.close()
