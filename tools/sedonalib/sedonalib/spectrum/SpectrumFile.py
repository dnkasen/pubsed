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

        self.read_data()


    # was hoping this would keep user from changing stuff...
    @property
    def mu(self):
        return self.mu


    def read_data(self):

        fin = h5py.File(self.filename,"r")
        self.x     = np.array(fin['nu'])
        self.L     = np.array(fin['Lnu'])
        self.click = np.array(fin['click'])

        self.t         = np.array(fin['time'])
        self.t_edges   = np.array(fin['time_edges'])
        self.mu        = np.array(fin['mu'])
        self.mu_edges  = np.array(fin['mu_edges'])
        self.phi       = np.array(fin['phi'])
        self.phi_edges = np.array(fin['phi_edges'])

        fin.close()

        self.n_x = len(self.x)
        self.n_times = len(self.t)
        self.n_mu = len(self.mu)
        self.n_phi = len(self.phi)

        print self.n_mu,self.n_phi
        # reshape so we always have 4 dimensions
        print self.L.shape
        if (self.n_mu == 1 and self.n_phi == 1):
            self.L = self.L[..., newaxis,newaxis]
        if (self.n_mu > 1 and self.n_phi == 1):
            self.L = self.L[..., newaxis]

        print self.L.shape


#        if (self.spec_units == 'angstrom'):
#            self.L  = self.L*self.x**2/pc.c/pc.cm_to_angs
#            self.x  = pc.c/self.x*pc.cm_to_angs

        if (self.time_units == 'day'):
            self.t       = self.t*pc.sec_to_day
            self.t_edges = self.t_edges*pc.sec_to_day


    def __str__(self):

        val = "Spectrum from file = " + self.filename
        #print ("n frequency = " + str(self.n_x))
        return val

    def switch_units(self,spec_units='angstroms'):

        new_spec_units = spec_units

    def get_bolometric_lc(self,angle_average=False,magnitudes=False):

        Lbol = np.zeros([self.n_times,self.n_mu,self.n_phi],dtype='d')
        Lave = np.zeros([self.n_times])
        print Lave.shape

        # integrate bolometric light curve
        for it,im,ip in np.ndindex(self.n_times,self.n_mu,self.n_phi):
            if (self.n_x == 1):
                Lbol[it,im,ip] = self.L[it,0,im,ip]
                Lave[it] += self.L[it,0,im,ip]
            else:
                Lbol[it,im,ip] = np.trapz(self.L[it,:,im,ip],x=self.x)
                Lave[it] += Lbol[it,im,ip]

        Lave = Lave/(1.0*self.n_mu*self.n_phi)

        if (magnitudes):
            Lbol = -2.5*np.log10(Lbol)+88.697425
            Lave = -2.5*np.log10(Lave)+88.697425

        if (angle_average):
            return self.t,Lave

        if (self.n_mu == 1 and self.n_phi == 1):
            return self.t,Lbol[:,0,0]
        if (self.n_mu > 1 and self.n_phi == 1):
            return self.t,Lbol[:,:,0]
        return self.t,Lbol

    def get_ABMag(self,band):
        """
        returns the AB magnitude light curve, for a given band
        """
        lum = np.zeros(self.time.size)
        dnu = np.ediff1d(self.nu,to_end=0)
        for i in range(self.time.size):
            sample_points = self.spec_tnu[i, :]/self.nu*(self.filter.transFunc_nu(band)(self.nu))
            lum[i] = inte.trapz(sample_points,x=self.nu)
            lum = lum/self.filter.getNormalization(band,self.nu) #gets Lnu(band)
            flx = lum/(4.*np.pi*(10.*3.0857e18)**2) #convert to flux at 10pc
            flx[np.where(flx==0.)] 	= 1e-99 #some small number so not taking log(0)
            mag = -2.5*np.log10(flx)-48.6 #get the AB magnitude
            mag[np.where(mag>0)] = 0. #set minimum mag to 0
            return mag
