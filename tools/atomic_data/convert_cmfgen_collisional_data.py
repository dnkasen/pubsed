import h5py
import numpy as np
import re
import matplotlib.pyplot as plt
import constants as pc




def add_col(out_fname,lev_fname,col_fname,species,ion,remove_J=True):

    print('\n-----------------')
    print(species + ' ' + ion)
    print('-----------------\n\n')

    fout = h5py.File(out_fname,'a')
    if (not str(species) in fout):
        species_group = fout.create_group(str(species))
    ion_group = fout.create_group(str(species) +'/' + str(ion))

    ##########################################
    ##### read level data
    fin = open(lev_fname,'r')
    # read off header
    line = fin.readline()
    while ('!Number of transitions' not in line):
        line = fin.readline()
    line = fin.readline()
    line = fin.readline()

    lev_c = []
    lev_E = []
    lev_g = []

    invcm_to_ev = 0.00012

    # read level data
    while (len(line.split()) > 1):
        data = line.split()
        config = data[0]
        if (remove_J):
            config = re.sub("[\[].*?[\]]", "", config)
        lev_c.append(config)
        lev_E.append(float(data[2])*invcm_to_ev)
        lev_g.append(float(data[1]))
        line = fin.readline()

    for i in range(len(lev_c)):
        print(lev_c[i],lev_g[i],lev_E[i])


    ######################################
    ## read collisional data
    fin = open(col_fname,'r')

    # read off header
    line = fin.readline()
    while ('Transition' not in line):
        line = fin.readline()

    # read temperature array (units = 10^4 K)
    data = line.split()
    data.pop(0)
    for i in range(len(data)):
        data[i] = data[i].replace('D','e')
        print(data[i])
    temps = np.array(data,dtype='d')

    # read data
    cnt = 0
    for line in fin:

        line = re.sub(' +-', '-', line)
        data = line.split()
        if (len(data) == 0): continue
        configs = data[0].split('-')
        if (remove_J):
            configs[0] = re.sub("[\[].*?[\]]", "", configs[0])
            configs[1] = re.sub("[\[].*?[\]]", "", configs[1])
        #configs[1] = configs[1].replace('-','')
        lev_l = lev_c.index(configs[0])
        lev_u = lev_c.index(configs[1])

        print(configs[0],configs[1],lev_l,lev_u,lev_g[lev_l])
        data.pop(0)
        coldata = np.array(data,dtype='d')
        coldata = 8.63e-08*coldata/temps**0.5/lev_g[lev_l]

        col_group = fout.create_group(str(species) +'/' + str(ion) + '/' + str(cnt))
        col_group.create_dataset("T", data=temps*1e4)
        col_group.create_dataset("C", data=coldata)
        col_group.create_dataset("lev_l",data=lev_l)
        col_group.create_dataset("lev_u",data=lev_u)

#        plt.plot(temps*1e4,coldata,'o')
#        T = temps*1e4
#        a = np.polyfit(T,coldata,4)
#        print(a)
#        afit = a[4] + a[3]*T + a[2]*T**2 + a[1]*T**3 + a[0]*T**4.0
#        plt.plot(T,afit)
#        plt.ion()
#        plt.yscale('log')
#        plt.show()
        print(species,ion,cnt)
#        j = input()
#        exit(0)
#        plt.clf()

        cnt += 1
    ion_group.create_dataset("n_data",data=cnt)

out_fname = 'cmfgen_col_data.h5'

base = '/Users/kasen//data/cmfgen_atomic_raw/'
col_fname = base + 'OXY/I/20sep11/oi_col'
lev_fname = base + 'OXY/I/20sep11/oi_osc_mchf'
species = '8'
ion     = '0'
add_col(out_fname,lev_fname,col_fname,species,ion)

col_fname = base + 'CA/II/30oct12/ca2col.dat'
lev_fname = base + 'CA/II/30oct12/ca2_osc_split.dat'
species = '20'
ion     = '1'
add_col(out_fname,lev_fname,col_fname,species,ion)

col_fname = base + 'FE/II/10sep16/fe2_col.dat'
lev_fname = base + 'FE/II/10sep16/fe2_osc'
species = '26'
ion     = '1'
add_col(out_fname,lev_fname,col_fname,species,ion,remove_J=False)
