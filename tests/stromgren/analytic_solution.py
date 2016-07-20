import pylab as py


kpc     = 3.08e21
pi      = 3.14159
clight  = 2.99792458e10    # speed of light (cm/s)
rydberg = 2.1798741e-11 # rydberg constant in ergs
k       = 1.380658e-16     # boltzmann constant (ergs/K)


rbox  = 6.6*kpc  
Ls    =  1.0933e38
ndens = 1e-3
temp  = 1e4


dr = 0.01*rbox
rr = py.arange(dr,2*rbox,dr)

# collisional ionization rate, approximately
zeta = rydberg/k/temp
f    = 2.7/zeta/zeta*temp**(-1.5)*py.exp(-1.0*zeta);
# recombination coefficient
alpha = 1.8e-13
# cross-section
sigma = 6.3e-18

tau = 0.0
xHI = 0*rr
for i in range(0,len(rr)):
 
    J = (Ls/(4.0*pi)**2/rr[i]**2)*py.exp(-tau)
    JN = 4*pi*J/rydberg*sigma/ndens

    b = (2*alpha + f + JN)/(alpha + f)
    c = 4*alpha/(alpha + f)
    det = b*b - c
    if (det < 0): det = 0
    xHI[i] = (b - py.sqrt(det))/2.0;


    tau += sigma*dr*ndens*xHI[i]

    print "{0:20.8e} {1:20.8e}".format(rr[i],xHI[i])

#py.plot(rr/rbox,xHI)
#py.yscale('log')
#py.ion()
#py.show()
#j = raw_input()
    
