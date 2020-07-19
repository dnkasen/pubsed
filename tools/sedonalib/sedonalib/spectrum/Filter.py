import os
import numpy as np
from scipy.interpolate import interp1d

class Filter:

    """
    Class to read and store filter transmisison functions
    """


    def __init__(self):

        self.Bands = {}
        self.FilterDatasets = [[]]

        datadir = os.environ['SEDONA_HOME'] + '/tools/sedonalib/sedonalib/data/filters/'
        with open(datadir + 'FILTER_LIST') as f:
            for line in f:
                dat = line.split()
                self.Bands[dat[1]] = int(dat[0])

        with open(datadir + 'allfilters.dat') as f:
            for line in f:
                if line.startswith('#'):
                    if self.FilterDatasets[-1] != ():
                        self.FilterDatasets.append([])
                else:
                    stripped_line = line.strip()
                    if stripped_line:
                        self.FilterDatasets[-1].append(stripped_line.split())

    # returns the raw filter transmission curve as {Wavelength (Angstrom), Transmission}
    def getTransmission(self,filt):
        filt_idx = self.Bands[filt]
        return np.array(self.FilterDatasets[filt_idx],dtype=np.float32)



    # returns an interpolation function for the transmission T(lambda)
    # from the raw data
    def transFunc(self,filt):
        wv,f = np.transpose(self.getTransmission(filt))
        interp = interp1d(wv,f,bounds_error=False,fill_value=0.0)
        return interp

    # returns an interpolation function for the transmission T(nu)
    def transFunc_nu(self,filt):
        wv,f = np.transpose(self.getTransmission(filt))
        nu = 3e10/(wv/1e8)
        interp = interp1d(nu,f,bounds_error=False,fill_value=0.0)
        return interp

    # get the filter boundaries
    def filtBounds(self,filt):
        return np.min(self.getTransmission(filt)[:,0]),np.max(self.getTransmission(filt)[:,0])

    # returns the normalization for the transmission
    def getNormalization(self,filt,lam):
        sample_points = self.transFunc(filt)(lam)/lam
        return np.trapz(sample_points,lam)

    # returns the normalization for the transmission
    def getNormalization_nu(self,filt,nu):
        sample_points = self.transFunc_nu(filt)(nu)/nu
        return np.trapz(sample_points,nu)

    def list(self):

        print('-----------------------------')
        print('List of Available Filters:')
        print('-----------------------------')
        for b in self.Bands:
            print(b)
        print('-----------------------------')
        print('See file FILTER_LIST for details/references')
        print('-----------------------------')
