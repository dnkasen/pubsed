import numpy as np

class SedonaBaseModel():


    @property
    def n_zones(self):
        return self.n_zones

    @property
    def n_elements(self):
        return len(self.elem_A)

    @property
    def dims(self):
        return self.dims

    @property
    def density(self):
        return self.dens

    @property
    def rho(self):
        return self.dens

    @property
    def ke(self):
        return self.kinetic_energy

    @property
    def composition(self):
        return self.comp

    def __init__(self):
        pass

    def __init__(self,name):
        pass

    def __str__(self):
        pass


    def set_constant_density(self,v):
        self.dens.fill(v)

    def set_constant_temperature(self,v):
        self.temp.fill(v)

    def set_constant_erad(self,v):
        self.erad.fill(v)

    def set_minimum_density(self,d):
        self.min_dens = d
        self.dens[self.dens < d] = d

    def set_minimum_temperature(self,t):
        self.min_temp = t
        self.temp[self.temp < t] = t

    def set_constant_composition(self,c,elements=None):

        if (elements is not None):
            self.set_elements(elements)

        nelem = len(self.elem_Z)
        if (len(c) != nelem):
            message = "compostion passed does not have same # of "
            message += "elements as model"
            raise ValueError(message)

        for i in range(nelem):
            if (len(self.dims) == 1):
                self.comp[:,i].fill(c[i])
            if (len(self.dims) == 2):
                self.comp[:,:,i].fill(c[i])
            if (len(self.dims) == 3):
                self.comp[:,:,:,i].fill(c[i])


    def set_density(self,d):

        if (np.isscalar(d)):
            self.set_constant_density(d)
            return

        if (d.shape != self.dims):
            message =  "dimensions of array not equal to "
            message += "grid dimensions = " + str(self.dims) + "\n"
            raise ValueError(message)

        self.dens = np.copy(d)


    def set_temperature(self,t):

        if (np.isscalar(t)):
            self.set_constant_temperature(t)
            return

        if (t.shape != self.dims):
            message =  "dimensions of array not equal to "
            message += "grid dimensions = " + str(self.dims) + "\n"
            raise ValueError(message)

        self.temp = np.copy(t)


    def set_erad(self,v):

        if (np.isscalar(v)):
            self.set_constant_erad(v)
            return

        if (v.shape != self.dims):
            message =  "dimensions of array not equal to "
            message += "grid dimensions = " + str(self.dims) + "\n"
            raise ValueError(message)

        self.erad = np.copy(v)


    def set_elements(self,elist):

        # wipe whatever elements might have been there
        self.elem_A = []
        self.elem_Z = []
        self.elem_list = []
        self.add_elements(elist)

    def add_elements(self,elist):

        if not isinstance(elist, (list, tuple)):
            elist = [elist]

        n_add = 0
        for e in elist:

            if (not isinstance(e, basestring)):
                e = str(e)

            if (e not in self.elem_list):
                self.elem_list.append(e)
                Z,A =  e.split(".")
                self.elem_Z.append(int(Z))
                self.elem_A.append(int(A))
                n_add += 1

        if (n_add != 0):
            nelem = len(self.elem_Z)

            if (self.comp is not None):
                self.comp.resize(self.dims + (nelem,))

    def add_element(self,e):
        self.add_elements(e)


    def check_base_model_validity(self):

        if (self.elem_A is None or self.elem_Z is None):
            print 'Error: no elements defined in composition'
            return False

        if (len(self.elem_Z) <= 0):
            print 'Error: no elements defined in composition'
            return False

        if (self.dens is None):
            print 'Error: density has not been set'
            return False

        if (self.temp is None):
            print 'Error: temperature has not been set'
            return False

        if (self.comp is None):
            print 'Error: composition has not been set'
            return False

        if (np.any(self.dens <= 0)):
            print 'Warning: zero or negative densities in model'

        if (np.any(self.temp <= 0)):
            print 'Warning: zero or negative temperatures in model'


        return True

    def write(self,outfile):
        pass


    def wait_for_key(self):
        j = raw_input("press any key to continue>")
