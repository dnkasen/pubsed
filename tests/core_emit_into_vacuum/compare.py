import os
import pylab as py
import h5py


def compare(pdf):
     
    ###########################################
    # compare the output
    ###########################################
    py.clf()


    h   = 6.6260755e-27    # planck's constant (ergs-s)
    c   = 2.99792458e10    # speed of light (cm/s)
    k   = 1.380658e-16     # boltzmann constant (ergs/K)
    sb  = 5.6704e-5        # stefan boltzman constant (ergs cm^-2 s^-1 K^-4) 
    pi  = 3.14159          # just pi
    py.ion()

    T   = 1e4
    L  = 1e43
    r0 = 0.5e15

    data = py.loadtxt('ray_00001')
    r = data[:,0] 
    trad = data[:,4]
    tgas = data[:,3]

    # dillution factor W
    rr = py.arange(r0,max(r),0.05*r0)
    w    = 0.5*(1 - (1 - r0**2/rr**2)**0.5)
    TW = (L*w/(4.0*pi*r0**2)/sb)**0.25

    py.plot(r,trad,'o',color='black')
    py.plot(r,tgas,'--',color='blue')
    py.plot(rr,TW,color='red',linewidth=2)
    py.legend(['sedona Trad','sedona Tgas','analytic solution'])
    py.title('core into vacuum test: radiation field')
    py.xlabel('radius (cm)')
    py.ylabel('radiation temperature (aT^4 = erad)')


    if (pdf != ''): pdf.savefig()
    else:
        py.show()
        j = raw_input()
        

    #------------------------------------------
    #compare output spectrum
    
#    py.clf()
    data = py.loadtxt('spectrum_1.dat')
    nu = data[:,0]
    y  = data[:,1]
    py.plot(nu,y,'o',color='black')
    # blackbody spectrum
    f = 2.0*h*nu**3/c**2/(py.exp(h*nu/k/T) - 1)
    f = f/(sb*T**4/pi)*L
    py.plot(nu,f,color='red',linewidth=2)

    py.legend(['sedona','analytic blackbody'])
    py.title('core into vacuum test: output spectrum')
    py.xlabel('frequency (Hz)')
    py.ylabel('Flux')
    py.yscale('log')
    py.xlim(0,3e15)
    py.ylim(1e25,1e29)

    if (pdf != ''): pdf.savefig()
    else:
        py.show()
        j = raw_input()


    # compare spectrum at 50
    py.clf()
    fin = h5py.File('grid_00001.h5','r')
    nu = py.array(fin['nu'])
    Jnu = py.array(fin['50/Jnu'])
    fin.close()
    py.plot(nu,Jnu,'o')

    # blackbody spectrum
    T = 1e4
    f = 2.0*h*nu**2.0/c**2/(py.exp(h*nu/k/T) - 1)*nu
    py.plot(nu,0.5*f,color='red',linewidth=2)
    
    py.legend(['sedona','analytic blackbody'])
    py.title('core into vacuum test: Jnu at zone=50')
    py.xlabel('frequency (Hz)')
    py.ylabel('Flux')
    py.yscale('log')
    py.xlim(0,3e15)
    py.ylim(5e-8,5e-4)

    
    if (pdf != ''): pdf.savefig()
    else:
        py.show()
        j = raw_input()
        
compare('')
