import h5py
import pylab as py




def add_element(z,fname,max_level,maxion):

    names = ['','h','he','li','be','b','c','n','o','f','ne','na','mg','al','si','p','s','cl','ar','k','ca','sc','ti','vn','cr','mn','fe','co','ni']

    if (z < maxion): maxion = z
 
    icm_to_ev = 1.23981e-4
    
    # open hdf5 file
    f = h5py.File(fname,'a')
    base = str(z) + "/"
    print base
    grp = f.create_group(base)
    
    # read ionization potential
    chi = []
    fin = open('ip/chianti.ip')
    for line in fin.readlines():
        if (line == ' -1\n'): break
        data = line.split()
        c = float(data[2])*icm_to_ev
        if (int(data[0]) == z): chi.append(c)
        if (len(chi) >= maxion): break
    fin.close()
    
    #---------------------------------------
    # read files
    #---------------------------------------
    E_lev = []
    i_lev = []
    c_lev = []
    g_lev = []
    ground = []
    line_l = []
    line_u = []
    line_A = []
    
    cnt = 0
   
    for j in range(1,maxion+1):


        ground.append(cnt)
        iname = names[z] + '_' + str(j)
        direc = names[z] + '/' + iname + '/'

        print direc + '.elvlc'
        #---------------------------------------
        # read level file
        #---------------------------------------
        try:
            fin = open(direc + iname + '.elvlc')
        except:
            # add in a single level
            E_lev.append(0)
            g_lev.append(1)
            c_lev.append('base')
            i_lev.append(j-1)
            cnt += 1
            continue
        
        lev_cnt = 0
        for line in fin:

            if (len(line.split()) < 2): break
            c_lev.append(line[7:37])

            # cut down to line
            print line
            line = line[52:len(line)]
            data = line.split()
            g_lev.append(float(data[0])*2 + 1)
            if (float(data[1]) >= 0):
                E_lev.append(float(data[1 ])*icm_to_ev)
                i_lev.append(j-1)
            else:
                E_lev.append(float(data[2])*icm_to_ev)
                i_lev.append(j-1)
            cnt += 1
            print float(data[1 ])*icm_to_ev, float(data[2 ])*icm_to_ev
            lev_cnt += 1
            if (lev_cnt >= max_level): break
                        
        print z,iname,lev_cnt
        fin.close()
        
        #---------------------------------------
        # read line file
        #---------------------------------------
        try:
            fin = open(direc + iname + '.wgfa')
        except:
            continue

        for line in fin:
            data = line.split()
            if (len(data) < 2): break
            ll = int(data[0]) - 1 + ground[j-1]
            lu = int(data[1]) - 1 + ground[j-1]

            if (int(data[1]) > max_level): continue

            line_l.append(ll)
            line_u.append(lu)
            line_A.append(float(data[4]))
            
    #    for i in range(len(E_lev)):
#        print i,c_lev[i],E_lev[i],g_lev[i],i_lev[i]
      

    print z,len(chi)

    # add extra ionization state
    E_lev.append(0)
    g_lev.append(1)
    i_lev.append(len(chi))
    chi.append(1000000)
    ground.append(len(E_lev)-1)


    # write out hdf5 file
    f[base].attrs["n_ions"]    = len(chi)
    f[base].attrs["n_levels" ] = int(len(E_lev))
    f[base].attrs["n_lines" ]  = len(line_A)

    f.create_dataset(base + "ion_chi",data = chi)
    f.create_dataset(base + "ion_ground",data = ground)
    f.create_dataset(base + "level_g",data = g_lev)
    f.create_dataset(base + "level_E",data = E_lev)
    f.create_dataset(base + "level_i",data = i_lev)
    f.create_dataset(base + "line_l",data = line_l)
    f.create_dataset(base + "line_u",data = line_u)
    f.create_dataset(base + "line_A",data = line_A)



outfile = 'chianti_atomdata_i4_l20.h5'


for i in range(1,28):
    add_element(i,outfile,800,4)
    




