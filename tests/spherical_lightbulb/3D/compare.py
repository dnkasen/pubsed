import os
import matplotlib.pyplot as py
import numpy as np
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

    # read data
    fin = h5py.File('plt_00001.h5')
    trad = np.array(fin['T_rad'])
    tgas = np.array(fin['T_gas'])
    r = np.array(fin['x'])
    #r = r - 0.5*(r[1] - r[0])
    nx = len(r)
    
    py.matshow(trad[:,:,nx/2])
    py.colorbar()
    py.title('T_rad (x-y plane)')
    if (pdf != ''): pdf.savefig()
    else:
        py.show()
        j = raw_input()
    py.clf()

    # dillution factor W
    rr = np.arange(r0,max(r),0.05*r0)
    w    = 0.5*(1 - (1 - r0**2/rr**2)**0.5)
    TW = (L*w/(4.0*pi*r0**2)/sb)**0.25

    py.plot(r,trad[:,nx/2,nx/2],'o',color='blue')
    py.plot(r,trad[nx/2,:,nx/2],'o',color='red')
    py.plot(r,trad[nx/2,nx/2,:],'o',color='black')

    #py.plot(r,tgas[nx/2,nx/2,:],'--',color='blue')
    py.ylim(4000,15000)
    py.plot(rr,TW,color='k',linewidth=2)
    py.plot(-1.0*rr,TW,color='red',linewidth=2)

    py.legend(['sedona Trad(x)','sedona Trad(y)','sedona Trad(z)','analytic solution'])
    py.title('core into vacuum test: radiation field')
    py.xlabel('radius (cm)')
    py.ylabel('radiation temperature (aT^4 = erad)')

    fin.close()

    if (pdf != ''): pdf.savefig()
    else:
        py.show()
        j = raw_input()
        

    #------------------------------------------
    #compare output spectrum
    
    py.clf()
    fin = h5py.File('spectrum_1.h5')
    nu   = np.array(fin['nu'],dtype='d')
    Lnu  = np.array(fin['Lnu'],dtype='d')
    py.plot(nu,Lnu[0,:],'o',color='black')
    # blackbody spectrum
    f = 2.0*h*nu**3/c**2/(np.exp(h*nu/k/T) - 1)
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

    #------------------------------------------
    # compare spectrum at zone 50
    py.clf()
    fin = h5py.File('grid_00001.h5','r')
    nu = py.array(fin['nu'])
    Jnu = py.array(fin['zone_50/Jnu'])
    rz = py.array(fin['x'])
    fin.close()
    py.plot(nu,Jnu,'o')

    # blackbody spectrum
    T = 1e4
    L  = 1e43
    f = 2.0*h*nu**2.0/c**2/(py.exp(h*nu/k/T) - 1)*nu
    f = L*f/(sb*T**4*4.0*pi*r0**2)
    W = 0.5*(1 - (1 - (r0/rz[50])**2)**0.5)
    py.plot(nu,W*f,color='red',linewidth=2)
    
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
        
if __name__=='__main__': compare('')
