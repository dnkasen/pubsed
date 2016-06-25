import numpy as np

# parameters of model
lam0  = 1215.0
vmax  = 2.0e9
vphot = 1.0e9
tau0  = 4.0
lam1  = 1000.0
lam2  = 1400.0
dlam  = 1.0
nx    = 200

dv = vmax/(nx*1.0)

c_light = 2.99e10
pi = 3.14159

for lam in np.arange(lam1,lam2,dlam):
    vz = (lam - lam0)/lam0*c_light
    flux = 0
    for vp in np.arange(0.0,vmax,dv):
        vr = (vz**2 + vp**2)**0.5

        # get optical depth
        tau = tau0
        # tau is zero in or behind photosphere
        if (vr > vmax):  tau = 0
        if (vr < vphot): tau = 0
        if (vp < vphot and vz > 0): tau = 0
        # get exponential of tau
        if (tau > 0): etau = np.exp(-tau)
        else: etau = 1

        # pure scattering source function (dilution function)
        if (vr < vphot or vr > vmax): S = 0
        else: S = 0.5*(1 - (1 - (vphot/vr)**2)**0.5)
        
        # photospheric intensity
        Iphot = 0
        if (vp < vphot): Iphot = 1*(lam0/lam)**2
      
        flux += 2*pi*vp*(Iphot*etau + S*(1 - etau))*dv

    norm = pi*vphot**2
    print lam,flux/norm
    
