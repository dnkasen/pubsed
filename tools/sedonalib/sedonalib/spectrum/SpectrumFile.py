import numpy as np
import h5py
import physical_constants as pc
from numpy import newaxis

class SpectrumFile():
    """
    Class for reading in and handling spectrum files output
    by the Sedona code

    Allowed units:
    spec_units: 'hz','angstrom','cm'

    The spectrum data is stored as a 4D array of the form
    L[n_times, n_freq, n_mu, n_phi]
    where L is the specific luminosity (ergs/s/Hz or ergs/s/Ang)

    """
    def __init__(self,name,spec_units='angstrom',time_units='day',verbose=True):

        self.verbose = verbose
        self.filename = name

        allowed_spec_units = ['angstrom','hz','cm','micron']
        self.default_spec_unit = 'angstrom'
        allowed_time_units = ['sec','day']

        self.spec_units = spec_units
        if (self.spec_units not in allowed_spec_units):
            if (verbose):
                print ('ERROR: unknown spectrum units; using default ')
            self.spec_units = 'angstrom'

        if (time_units == "s"):
            time_units = "sec"
        if (time_units == "days"):
            time_units = "day"
        self.time_units = time_units
        if (self.time_units not in allowed_time_units):
            if (verbose):
                print ('ERROR: unknown time units; using default = day')
            self.time_units = 'day'

        self.Lbol = None
        self.read_data()


    # was hoping this would keep user from changing stuff...
    @property
    def mu(self):
        return self.mu

    @property
    def shape(self):
        return self.L.shape

    def read_data(self):

        fin = h5py.File(self.filename,"r")
        self.x     = np.array(fin['nu'])
        self.L     = np.array(fin['Lnu'])
 #       self.click = np.array(fin['click'])

        self.t         = np.array(fin['time'])
        try:
            self.mu  = np.array(fin['mu'])
        except:
            self.mu  = np.zeros(1)
        try:
            self.phi       = np.array(fin['phi'])
        except:
            self.phi = np.zeros(1)

    #    self.phi_edges = np.array(fin['phi_edges'])
   #     self.mu_edges  = np.array(fin['mu_edges'])
  #      self.t_edges   = np.array(fin['time_edges'])

        fin.close()

        self.n_x = len(self.x)
        self.n_times = len(self.t)
        self.n_mu = len(self.mu)
        self.n_phi = len(self.phi)
#        self.n_phi = 1

        # reshape so we always have 4 dimensions
        if (self.n_mu == 1 and self.n_phi == 1):
            self.L = self.L[..., newaxis,newaxis]
        if (self.n_mu > 1 and self.n_phi == 1):
            self.L = self.L[..., newaxis]


        if (self.spec_units == 'angstrom'):
            for it,im,ip in np.ndindex(self.n_times,self.n_mu,self.n_phi):
                # change units of luminosity to erg/s/A
                newL = self.L[it,:,im,ip]*self.x**2/pc.c/pc.cm_to_angs
                # reverse the order fo the array
                self.L[it,:,im,ip]  = newL[::-1]

            # change x units to wavelength and flip order
            newx  = pc.c/self.x*pc.cm_to_angs
            self.x = newx[::-1]

        if (self.time_units == 'day'):
            self.t       = self.t*pc.sec_to_day
 #           self.t_edges = self.t_edges*pc.sec_to_day


    def __str__(self):

        val =  "Spectrum from file = " + self.filename
        val += "\n\ntime pts = " + str(len(self.t))
        val += "\n freq pts = " + str(len(self.x))
        val += "\n   mu pts = " + str(len(self.mu))
        val += "\n  phi pts = " + str(len(self.phi))

        return val

    def switch_units(self,spec_units='angstroms'):

        new_spec_units = spec_units

    def get_band_lc(self,wrange,view=None,angle_average=False,magnitudes=False):

        Lband = np.zeros([self.n_times,self.n_mu,self.n_phi],dtype='d')

        b = (self.x >= wrange[0])*(self.x <= wrange[1])
        if (sum(b) == 0):
            message = "Wavelength range for band light curve not in spectrum range"
            raise ValueError(message)

        # integrate bolometric light curve
        for it,im,ip in np.ndindex(self.n_times,self.n_mu,self.n_phi):
            Lband[it,im,ip] = np.trapz(self.L[it,b,im,ip],x=self.x[b])

        Lband = Lband/(wrange[1] - wrange[0])

        if (self.n_mu == 1 and self.n_phi == 1):
            return self.t,Lband[:,0,0]
        if (self.n_mu > 1 and self.n_phi == 1):
            return self.t,Lband[:,:,0]
        if (self.n_mu == 1 and self.n_phi > 1):
            return self.t,Lband[:,0,:]
        return self.t,Lband

    def get_bolometric_lc(self,view=None,angle_average=False,magnitudes=False):

        # compute bolometric light curve
        if (self.Lbol is None):
            self.Lbol = np.zeros([self.n_times,self.n_mu,self.n_phi],dtype='d')
            self.Lave = np.zeros([self.n_times])

        # integrate bolometric light curve
        for it,im,ip in np.ndindex(self.n_times,self.n_mu,self.n_phi):
            if (self.n_x == 1):
                self.Lbol[it,im,ip] = self.L[it,0,im,ip]
                self.Lave[it] += self.L[it,0,im,ip]
            else:
                self.Lbol[it,im,ip] = np.trapz(self.L[it,:,im,ip],x=self.x)
                self.Lave[it] += self.Lbol[it,im,ip]

                self.Lave = self.Lave/(1.0*self.n_mu*self.n_phi)

        Lbol = self.Lbol
        Lave = self.Lave

        thisL = self.Lbol
        if (angle_average):
            thisL = self.Lave

        if (magnitudes):
            thisL = -2.5*np.log10(thisL)+88.697425

        if (self.n_mu == 1 and self.n_phi == 1):
            return self.t,thisL[:,0,0]
        if (self.n_mu > 1 and self.n_phi == 1):
            return self.t,thisL[:,:,0]

        # interpolate to particular viewing angle
        if (view is not None):
            mu  = view[0]
            phi = view[1]
            import bisect

            imu = bisect.bisect_left(self.mu,mu)
            i1 = imu-1
            i2 = imu
            if (i2 == len(self.mu)): i2 = imu-1
            m1 = self.mu[i1]
            m2 = self.mu[i2]
            dm = mu - m1
    #        print i1,i2,m1,m2,mu


            jphi = bisect.bisect_left(self.phi,phi)
            nphi = len(self.phi)
            if (jphi == 0):
                j1 = nphi-1
                j2 = 0
                p1 = self.phi[j1] - 2.0*np.pi
                p2 = self.phi[j2]
            elif (jphi == nphi):
                j1 = nphi-1
                j2 = 0
                p1 = self.phi[j1]
                p2 = self.phi[j2] + 2.0*np.pi
            else:
                j1 = jphi -1
                j2 = jphi
                p1 = self.phi[j1]
                p2 = self.phi[j2]

            print p1,phi,p2, jphi,j1,j2

            p1 = self.phi[j1]
            p2 = self.phi[j2]
            if (jphi == len(self.phi)):
                j2 = jphi - 1
                j1 = 0
                p1 = np.pi*2.0 + self.phi[0]


            dp = phi - p1

            Llo = self.Lbol[:,i1,j1]
            Lhi = self.Lbol[:,i2,j1]
            if (m2 == m1):
                L1 = Llo
            else:
                L1  =   Llo + (Lhi - Llo)*dm/(m2 - m1)

            Llo = self.Lbol[:,i1,j2]
            Lhi = self.Lbol[:,i2,j2]
            if (m2 == m1):
                L2 = Llo
            else:
                L2  = Llo + (Lhi - Llo)*dm/(m2 - m1)

            if (p2 == p1):
                Lint = L1
            else:
                Lint = L1 + (L2 - L1)*dp/(p2 - p1)

            thisL = Lint


        return self.t,thisL




    def get_spectrum(self,time=None,mu=None,phi=None,interpolate=True,angle_average=False):

        import bisect

        nt = len(self.t)
        nmu  = len(self.mu)
        nphi = len(self.phi)
        nx   = len(self.x)

        # for just a snapshot 1D spectrum
        if (nt == 1 and nmu == 1 and nphi == 1):
            return self.x, self.L[0,:,0,0]

        smp = np.zeros((nx,nmu,nphi))

        if (interpolate and nt > 1):

            # interpolate in time
            indt = bisect.bisect_left(self.t,time)
            i1 = indt-1
            i2 = indt
            t1 = self.t[i1]
            t2 = self.t[i2]
            dt = time - t1
            for i in range(nmu):
                for j in range (nphi):
                    s1 = self.L[i1,:,i,j]
                    s2 = self.L[i2,:,i,j]
                    smp[:,i,j] = s1 + (s2 - s1)*dt/(t2 - t1)

        else:
            if (nt == 1):
                smp = self.L[0,:,:,:]
            else:
                indt = bisect.bisect(self.t,time)
                smp = self.L[indt,:,:,:]

            # interpolate in mu, if wanted
#            if (mu is not None):
#                imu = bisect.bisect_left(self.mu,mu)
#                i1 = imu-1
#                i2 = imu
#                m1 = self.mu[i1]
#                m2 = self.mu[i2]
#                dm = mu - m1
#                for j in range (nphi):
#                    s1 = smp[:,i1,j]
#                    s2 = smp[:,i2,j]
#                    sm[:,j] = s1 + (s2 - s1)*dm/(m2 - m1)

        if (angle_average):
            Fave = np.zeros(self.n_x)
            for im,ip in np.ndindex(self.n_mu,self.n_phi):
                Fave += smp[:,im,ip]
            Fave /= 1.0/(1.0*self.n_mu*self.n_phi)
            return self.x, Fave

        if (nmu == 1 and nphi == 1):
            return self.x,smp[:,0,0]
        elif (nmu == 1):
            return self.x,smp[:,0,:]
        else:
            return self.x,smp[:,:,:]


    def get_nuLnu_lc(self,wrange,magnitudes=False):
        b = (self.x >= wrange[0])*(self.x <= wrange[1])
        new_x = self.x[b]
        lc = np.zeros(self.t.size)
        for i in range(self.t.size):
            lc[i] = np.trapz(self.L[i,b,0,0],x=new_x)


        if (magnitudes):
            lc = lc/(wrange[1] - wrange[0])
            
        else:
            lc = lc*0.5*(wrange[0] + wrange[1])/(wrange[1] - wrange[0])

        return self.t,lc


    def get_ABMag(self,band):
        """
        returns the AB magnitude light curve, for a given band
        """
        lum = np.zeros(self.t.size)
        dnu = np.ediff1d(self.nu,to_end=0)
        for i in range(self.time.size):
            sample_points = self.spec_tnu[i, :]/self.nu*(self.filter.transFunc_nu(band)(self.nu))
            lum[i] = inte.trapz(sample_points,x=self.nu)
            lum = lum/self.filter.getNormalization(band,self.nu) #gets Lnu(band)
            flx = lum/(4.*np.pi*(10.*3.0857e18)**2) #convert to flux at 10pc
            flx[np.where(flx==0.)] 	= 1e-99 #some small number so not taking log(0)
            mag = -2.5*np.log10(flx)-48.6 #get the AB magnitude
            mag[np.where(mag>0)] = 0. #set minimum mag to 0
            return self.t,mag
