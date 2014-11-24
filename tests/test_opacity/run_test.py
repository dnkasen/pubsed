import os
import pylab as py

def run(pdf,plotup):

    ###########################################
    # run the code
    ###########################################
    os.system("./test_opacity")
       
    ###########################################
    # compare the output
    ###########################################
    py.clf()

    temp,ion = py.loadtxt("H_ionization.txt",unpack=1)

    # constants
    h   = 6.6260755e-27    # planck's constant (ergs-s)
    k   = 1.380658e-16     # boltzmann constant (ergs/K)
    k_ev = 8.6173324e-5    # boltzmann constant (ev/K)
    pi  = 3.14159          # just pi
    m_e = 9.10938188e-28   # mass of electron (g)
    m_p = 1.67262158e-24   # mass of proton (g)


    # analytic ionization
    ndens = 1e-13/m_p
    lam = (h**2/(2*pi*m_e*k*temp))**0.5
    phi = 2/lam**3/2.0*py.exp(-13.605/k_ev/temp)
    b = 1.0*phi/ndens + 2
    f = -1*(-1.0*b + (b**2 - 4)**0.5)/2.0
    py.plot(temp,ion,'o',color='black',markersize=8,markeredgewidth=2,markerfacecolor='none')
    py.plot(temp,1-f,color='black')
    py.ylim(-0.1,1.1)
    py.xlabel('temperature (K)')
    py.ylabel('ionization fraction')
    py.title('hydrogen ionization test')
    py.legend(['sedona opacity module','analytic'],loc=4)

    pdf.savefig()
    if (plotup):
        py.ion()
        py.show()
        j = raw_input()

    ##############################
    # Iron ionization
    py.clf()

    temp,ion = py.loadtxt("Fe_ionization.txt",unpack=1)

    py.ylim(-0.1,10.1)
    py.xlabel('temperature (K)')
    py.ylabel('ionization fraction')
    py.title('iron ionization test')
    py.legend(['sedona opacity module','analytic'],loc=4)
    py.plot(temp,ion,'o',color='black',markersize=8,markeredgewidth=2,markerfacecolor='none')
    pdf.savefig()
    if (plotup):
        py.ion()
        py.show()
        j = raw_input()

  ##############################                                                                          
  # Iron ionization                                                                                       
    py.clf()
