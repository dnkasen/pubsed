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
    
    # compare spectrum 
    py.clf()
    fin = h5py.File('grid_00001.h5','r')
    nu = py.array(fin['nu'])
    Jnu = py.array(fin['zone_0/Jnu'])
    opac = py.array(fin['zone_0/opacity'])
    emis = py.array(fin['zone_0/emissivity'])
    Snu = py.nan_to_num(emis/opac)
    fin.close()
    py.plot(nu,Jnu,'o')

    # blackbody spectrum
    T = 1e5
    f = 2.0*h*nu**2.0/c**2/(py.exp(h*nu/k/T) - 1)*nu
    py.plot(nu,f,color='red',linewidth=2)
    
    py.legend(['sedona','analytic blackbody'])
    py.title('one zone eq test: radiation field')
    py.xlabel('frequency (Hz)')
    py.ylabel('Flux')
    py.yscale('log')
    py.xscale('log')
    py.xlim(1e14,5e16)
    py.ylim(1e-4,5)

    
    if (pdf != ''): pdf.savefig()
    else:
        py.show()
        j = raw_input()



    py.clf()

    fin = h5py.File('levels_00001.h5','r')
    b = py.array(fin['zone_0/Z_1/level_departure'])
    py.plot(b,'o',color='black')
    py.plot([0,len(b)+1],[1,1],color='black')
    py.title('one zone eq test: NLTE departure coefficients')
    py.xlabel('level #',size=15)
    py.ylabel('departure coefficient',size=15)

    for i in range(len(b)):
        if (b[i] > 1.2 or b[i] < 0.8):
            py.plot([i,i],[-1,10],color='red')

    py.ylim(0.8,1.2)


    fin.close()
    if (pdf != ''): pdf.savefig()
    else:
        py.show()
        j = raw_input()




if __name__=='__main__': compare('')
