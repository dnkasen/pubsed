import string
import numpy  as np
import h5py
import matplotlib.pyplot as plt
import bisect

def myshow():
    plt.ion()
    plt.show()
    j = raw_input('press any key to quit>')


def get_plt_time(fname):
    fin = h5py.File(fname,'r')
    t = np.array(fin['time'])
    return t[0]/3600./24.

def get_flux_mean_opacity(fname,z):
    fin  = h5py.File(fname,'r')
    Jnu  = np.array(fin['zonedata/' + str(z) + '/Jnu'])
    opac = np.array(fin['zonedata/' + str(z) + '/opacity'])
    nu   = np.array(fin['nu'])
    mean = np.trapz(opac*Jnu,nu)
    norm = np.trapz(Jnu,nu)
    fin.close()
    return mean/norm

def get_flux_mean_opacity(fname):

    fin  = h5py.File(fname,'r')
    r    = np.array(fin['r'])
    mean = 0*r
    for z in range(len(r)):
        Jnu  = np.array(fin['zonedata/' + str(z) + '/Jnu'])
        opac = np.array(fin['zonedata/' + str(z) + '/opacity'])
        nu   = np.array(fin['nu'])
        Jop  = np.trapz(opac*Jnu,nu)
        norm = np.trapz(Jnu,nu)
        if (norm == 0): 
            mean[z] = 0
        else:
            mean[z] = Jop/norm
    fin.close()
    return r,mean

def get_flux_mean_tau(fname):
  
    r,k = get_flux_mean_opacity(fname)
    fin  = h5py.File(fname,'r')
    rho = np.array(fin['rho'])
    alpha = k*rho
    fin.close()
    tau = 0*r
    for i in range(len(r)-2,-1,-1):
        dr = r[i+1] - r[i]
        a  = 0.5*(alpha[i+1] + alpha[i])
        tau[i] = tau[i+1] + dr*a
    return r,tau

def get_blackbody_lam(lam,T):

    h   = 6.6260755e-27    # planck's constant (ergs-s)
    c   = 2.99792458e10    # speed of light (cm/s)
    k   = 1.380658e-16     # boltzmann constant (ergs/K)
    pi  = 3.14159          # just pi

    zeta = h*c/(lam*1e-8)/k/T
    B = 2*h*c*c/lam**5/(np.exp(zeta)-1)

    # make per angstroms
    return B*1e-8


def get_blackbody_nu(nu,T):


    h   = 6.6260755e-27    # planck's constant (ergs-s)
    c   = 2.99792458e10    # speed of light (cm/s)
    k   = 1.380658e-16     # boltzmann constant (ergs/K)

    zeta = h*nu/k/T
    B = 2*h*nu**3/c**2.0/(np.exp(zeta)-1)
    return B


def get_bolometric(fname):

    fin = h5py.File(fname,'r')
    time = np.array(fin['time'],dtype='d')
    Lnu  = np.array(fin['Lnu'],dtype='d')
    nu   = np.array(fin['nu'],dtype='d')
    bol = 0*time
    for i in range(0,len(bol)):
        bol[i] = np.trapz(Lnu[i,:],nu)
    fin.close()
    time = time/3600.0/24.0
    return time,bol


def get_nuFnu_band(fname,band):


    nu1 = 3e10*1e8/band[1]
    nu2 = 3e10*1e8/band[0]
    dnu = nu2 - nu1
    print nu1,nu2

    fin = h5py.File(fname,'r')
    time = np.array(fin['time'],dtype='d')
    Lnu  = np.array(fin['Lnu'],dtype='d')
    nu   = np.array(fin['nu'],dtype='d')
    lc = 0*time
    nuc = nu[(nu >= nu1)*(nu <= nu2)]
    print nu
    for i in range(0,len(lc)):
        f = Lnu[i,:]
        f = f[(nu >= nu1)*(nu <= nu2)]
        lc[i] = np.trapz(f,nuc)
    fin.close()
    time = time/3600.0/24.0
    lc = lc*0.5*(nu2 + nu1)/dnu
    return time,lc


def get_spectrum(fname,t0):

    fin = h5py.File(fname,'r')
    time = np.array(fin['time'],dtype='d')
    Lnu  = np.array(fin['Lnu'],dtype='d')
    nu   = np.array(fin['nu'],dtype='d')
    time = time/3600.0/24.0

    indt = bisect.bisect(time,float(t0))
    print 'time = ',time[indt]

    # convert to Llam
    Llam = Lnu*nu**2/2.99e10/1e8
    lam  = 2.99e10/nu*1e8

    return lam,Llam[indt,:]


def timesplit_specfile(fname,imu):
    flux,time,wave,mu = read_specfile(fname)
    for i in range(len(time)):
        day = "{0:.2f}".format(time[i])
        if (time[i] < 10): day = "0{0:.2f}".format(time[i])
        name = 'spec_d' + day + '_m' + str(imu) + '.dat'
        if (imu < 10): 
            name = 'spec_d' + day + '_m0' + str(imu) + '.dat'

        fin = open(name,'w')
        for j in range(len(wave)):
            fin.write(str(wave[j]) + ' ' + str(flux[i,j,imu]) + '\n')
        fin.close()



def read_specfile(fname):

    try: f = open(fname,'r')
    except: 
        print "can't open file: " + fname
        exit()

    # get dimensions 
    line = f.readline()
    terms = string.split(line)
    nt = int(terms[0])
    nw = int(terms[1])
    nm = int(terms[2])
    np = int(terms[3])

    time = numpy.ndarray(nt)
    wave = numpy.ndarray(nw)
    mu   = numpy.ndarray(nm)
    phi  = numpy.ndarray(np)
    flux = numpy.ndarray((nt,nw,nm))
    bol      = numpy.zeros((nt,nm))
    bol_ave  = numpy.zeros(nt)
    ave_spec = numpy.zeros((nt,nw))

    f_ind = 2
    if (nm > 1): f_ind = f_ind + 1
    if (np > 1): f_ind = f_ind + 1


    # read in the data
    for im in range(nm):
        for it in range(nt):
            for iw in range(nw):
                line = f.readline()
                terms = string.split(line)
                time[it] = float(terms[0])/3600./24.
                wave[iw] = float(terms[1])
                mu[im]   = float(terms[2])
                flux[it,iw,im] = float(terms[f_ind])

    return (flux,time,wave,mu)


