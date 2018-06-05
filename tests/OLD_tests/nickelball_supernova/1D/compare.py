import os
import pylab as py
import h5py
import bisect
import smooth

def compare(pdf):
     
    ###########################################
    # compare the output
    ###########################################
    py.clf()

    fin = h5py.File('spectrum.h5','r')
    nu    = py.array(fin['nu'],dtype='d')
    times = py.array(fin['time'])
    Lnu   = py.array(fin['Lnu'],dtype='d')
    # calculate bolometric LC
    nt  = len(times)
    nnu = len(nu)
    bol = py.zeros(nt,dtype='d')
    for it in range(nt):
        bol[it] = py.trapz(Lnu[it,:],x=nu)
#        for j in range(1,len(nu)):
 #           bol[it] += Lnu[it,j]*(nu[j] - nu[j-1]) #py.trapz(Lnu[it,:],x=nu)

#    py.plot(nu,Lnu[25,:])
#    py.show()
    py.plot(times/3600.0/24.0,bol,linewidth=4,color='k')

    xc,rc,lc = py.loadtxt('../comparefiles/sedonabox_nickelball_lc.dat',usecols=[0,1,2],skiprows=1,unpack=1)
    py.plot(xc/3600.0/24.0,lc,'o',color='black')
    py.plot(xc/3600.0/24.0,rc)
  #  py.yscale('log')
    py.ylim(4e42,6e43)
    py.xlim(0,50)

    # gamma-ray deposition
    tdep = []
    gdep = []
    rdep = []
    for j in range(1,94):
        ray = 'ray_00'
        if (j< 100): ray = ray + '0' 
        if (j < 10): ray = ray + '0' 
        ray = ray + str(j) + '.dat'
        fin = open(ray,'r')
        line = fin.readline()
        tdep.append(float(line.split()[3])/3600.0/24.0)
        sum = 0
        rsum = 0
        r0 = 0 
        line = fin.readline()
        for line in fin:
            data = line.split()
            r1 = float(data[0])
            vol = 4.0*3.14159/3.0*(r1**3 - r0**3)
            sum += float(data[5])*vol
            rsum += float(data[6])*vol
            r0 = r1
        gdep.append(sum)
        rdep.append(rsum)
    py.plot(tdep,gdep,'o')
    py.plot(tdep,rdep,'o')



    if (pdf != ''): pdf.savefig()
    else:
        py.ion()
        py.show()
        j = raw_input()
        
    #------------------------------------------
    #compare  spectrum at day 20

    py.clf()

    # time of spectrum we want to plot
    tshow = 20.*24*3600.

    # read data from hdf5 file
    f = h5py.File('../comparefiles/sedonabox_nickelball.h5', 'r')
    strt_time = f['strt_time'][()]
    stop_time = f['stop_time'][()]
    nu_grid   = f['nu_grid'  ][:]
    optical   = f['optical'  ][:,0,0,:]
    f.close()
    # number of time points and frequency points
    nt    = optical.shape[0]
    nfreq = optical.shape[1]
    # determine time spacing and index of time slice
    dt = (stop_time-strt_time) / float(nt)
    it = int((tshow-strt_time)/dt)
    # get center of frequency bins
    nuc = 0.5*(nu_grid[1:] + nu_grid[0:-1])
    # get width of frequency bins
    dnu = nu_grid[1:] - nu_grid[0:-1]
    # get energy in each freq. bin at this time (units ergs)
    Lbox = optical[it,:]
    # convert to specific luminosity (ergs/sec/Hz)
    Lnu_box = Lbox/dnu/dt

    # convert frequency to wavelength
    c = 2.99e10
    lam = c/nuc*1e8
    Llam = Lnu_box*nuc/lam
    # plot up the result
    py.plot(lam, Llam,linewidth=3)

    py.xlim(100,1e4)
    py.xlabel('wavelength (Angstroms)')
    py.ylabel('specific luminosity (ergs/sec/Angstrom)')
    py.title('t = %.1f days' %(tshow/3600.0/24.0))

    indt = bisect.bisect(times,tshow)
    Lnu = Lnu*nu**2/2.99e10/1e8
    nu  = 2.99e10/nu*1e8

    py.plot(nu,smooth.smooth(Lnu[indt,:],1))

    if (pdf != ''): pdf.savefig()
    else:
        py.show()
        j = raw_input()


if __name__=='__main__': compare('')
