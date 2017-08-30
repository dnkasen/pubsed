#!/usr/bin/env python
 
import bisect
import optparse
import numpy as np
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

    T   = 1.0e4
    L  = 1e43
    r0 = 0.5e15



    fin = h5py.File('zone_00001.h5','r')
    dr   = py.array(fin['dr'])
    z    = py.array(fin['z'])
    p    = py.array(fin['r'])
    Trad = py.array(fin['T_rad'])

    
    py.matshow(Trad)    
    py.title('2D core into vacuum test: radiation field')
    
    if (pdf != ''): pdf.savefig()
    else:
    	py.show()
    	j = raw_input()
    py.clf()


    for j in (25,50,100):
	    py.plot(p,Trad[:,j],'o')

	    # dillution factor W
	    rr = (p*p + z[j]*z[j])**0.5
	    w    = 0.5*(1 - (1 - r0**2/rr**2)**0.5)
	    TW = (L*w/(4.0*pi*r0**2)/sb)**0.25
	    py.plot(p,TW,color='k',linewidth=3)
    
    py.legend(['sedona Trad','analytic solution'])
    py.title('2D core into vacuum test: radiation field')
    py.xlabel('radius (cm)')
    py.ylabel('radiation temperature (aT^4 = erad)')

    if (pdf != ''): pdf.savefig()
    else:
    	py.show()
    	j = raw_input()
    fin.close()

    #------------------------------------------
    #compare output spectrum
    
    py.clf()
    fin = h5py.File('spectrum_1.h5','r')
    nu  = py.array(fin['nu'])
    Lnu = py.array(fin['Lnu'])
    nmu = (Lnu.shape)[2]
    for i in range(nmu):
	    py.plot(nu,Lnu[0,:,i],'o') #,color=cmap[i])
    # blackbody spectrum
    f = 2.0*nu*h*nu**2.0/c**2/(py.exp(h*nu/k/T) - 1)
    f = f/(sb*T**4/pi)*L
    py.plot(nu,f,color='black',linewidth=5)
    py.plot(nu,f,color='white',linewidth=3)


    py.legend(['sedona','analytic blackbody'])
    py.title('2D core into vacuum test: output spectrum')
    py.xlabel('frequency (Hz)')
    py.ylabel('Flux')
    py.yscale('log')
    py.xscale('log')
    py.xlim(2e13,4e15)
    py.ylim(1e24,3e28)

    if (pdf != ''): pdf.savefig()
    else:
        py.show()
        j = raw_input()

    #------------------------------------------


if __name__=='__main__': compare('')
