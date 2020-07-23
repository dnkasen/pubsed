import numpy as np
import h5py

import warnings
import inspect

def SedonaPltFile(name):
    """
    Factory function for loading Sedona plt files. See help(SedonaPltFileBase) or call
    help on the object itself for a description of the available methods.
    """

    data, dim = load_plt_data(name)

    if dim == 1:
        return Sedona1DSpherePltFile(name, data=data)
    elif dim == 2:
        return Sedona2DCylnPltFile(name, data=data)
    elif dim == 3:
        return Sedona3DCartPltFile(name, data=data)
    else:
        warnings.warn("Unrecognized dimensionality: " + str(dim) +
                ". Defaulting to generic SedonaPltFile base class.")
        return SedonaPltFileBase(name, data)

def load_plt_data(name):
    """
    Load data from a plt file, and return the resulting data dictionary and the
    dimensionality of the dataset (as a tuple, in that order).
    """

    ds      = h5py.File(name, 'r')
    keys    = ds.keys()
    data    = {key: np.array(ds[key], dtype = np.float64) for key in keys}
    ds.close()

    dim = len(data['rho'].shape)

    print("Loaded {}. Found it to be a {}-dimensional setup.".format(repr(name), dim))

    return (data, dim)

def derived_field(func):
    """ Decorator for flagging a method as a derived field. Will also turn it into a property. """

    func._derived_field = True
    return property(func)

def fetch_derived_fields(obj):
    """ Returns a generator of all of the derived fields the argument possesses. """

    for attr_name in dir(type(obj)):
        
        attr = inspect.getattr_static(obj, attr_name)
        if not isinstance(attr, property):
            continue
        if getattr(attr.fget, '_derived_field', False):
            yield attr_name

class SedonaPltFileBase():

    def __init__(self, name, autoload=True, data=None):

        self.name = name

        if data is not None:
            # If data was supplied, do not open file
            self.data = data
            self.dim = len(data['rho'].shape)
            self.data_status = "loaded"
            self.data_keys = self.data.keys()
        elif autoload:
            # Load data from file
            self.load_data()
        else:
            # Delay loading until later
            self.data = self.dim = self.data_keys = None
            self.data_status = "unloaded"

        self.derived_fields = set(fetch_derived_fields(self))
        self.fields = set(self.data_keys) | self.derived_fields

        #
        # What follows should be set by the extending subclass
        #

        # Tuple of coordinates (e.g. ('r', 'z'))
        self.coords    = ()
        # Velocity data keys, cleared upon loading data
        self._vel_keys = []

        # Default plot variables
        self.default_x_key = None
        self.default_y_key = None
        self.default_z_key = None

        if type(self) == SedonaPltFileBase:
            warnings.warn("Not setting any coordinates or default options for the plot")

    def __getitem__(self, key):
        """
        Dictionary-like data access. Also provides access to derived fields,
        while using something of the form 'self.data[key]' does not.
        """

        assert self.data is not None, "Must call load_data() before attempting data access"

        if key in self.data:
            return self.data[key]
        elif key in self.derived_fields:
            return getattr(self, key)
        else:
            raise KeyError("Field {} is not in dataset {}".format(repr(key), repr(self.name)))

    def load_data(self):
        """ Load the data contained within the plt file. """

        self.data, self.dim = load_plt_data(name)
        self.data_status = "loaded"
        self.data_keys = self.data.keys()
        self._vel_keys = []

    def vel_keys(self):
        """ A list of the velocities in the dataset, in (r, x, y, z) order. """
        
        if not self._vel_keys:

            order = ('r', 'x', 'y', 'z')
            vel_keys = filter(lambda k: k.startswith('vel'), self.data_keys)
            vel_keys = filter(lambda k: k[3:] in order, vel_keys)
            vel_keys = sorted(vel_keys, key=lambda k: order.index(k[3:]))
            self._vel_keys = vel_keys

        return self._vel_keys

    ###############################################
    # Derived fields
    ###############################################

    @derived_field
    def cell_vol(self):
        """ Array containing the volume of each cell. """

        raise NotImplementedError

    @derived_field
    def cell_mass(self):
        """ Array containing the mass of the material in each cell. """

        return self.cell_vol * self['rho']

    @derived_field
    def vtot(self):
        """ Array containing the total velocity in each cell. """

        vkeys = self.vel_keys()
        offset = max(len(vkeys) - self.dim, 0)
        if len(vkeys) != self.dim:
            warnings.warn("Using {} to calculate total velocity.".format(vkeys[offset:]))
        return sum(self[vi]**2 for vi in vkeys[offset:])**0.5

    @derived_field
    def kinetic_energy(self):
        """ Array containing the translational kinetic energy in each cell. """

        vkeys = self.vel_keys()
        offset = max(len(vkeys) - self.dim, 0)
        if len(vkeys) != self.dim:
            warnings.warn("Using {} to calculate kinetic energy.".format(vkeys[offset:]))
        return 0.5 * self.cell_mass * sum(self[vi]**2 for vi in vkeys[offset:])

    @derived_field
    def total_mass(self):
        """ The total mass contained within the domain. """

        return np.sum(self.cell_mass)

class Sedona1DSpherePltFile(SedonaPltFileBase):

    def __init__(self, name, autoload=True, data=None):

        super().__init__(name, autoload, data)

        self.default_x_key = 'r'
        self.default_y_key = 'rho'
        self.coords = ('r',)

    def load_data(self):

        super().load_data()

        assert self.dim == 1, "Attempting to load 1D spherical data that is not 1D"

    @derived_field
    def cell_vol(self):

        dr, = self['dr']
        return 4*np.pi/3 * (self['r']**3 - (self['r'] - dr)**3)

class Sedona2DCylnPltFile(SedonaPltFileBase):

    def __init__(self, name, autoload=True, data=None):

        super().__init__(name, autoload, data)

        self.default_x_key = 'r'
        self.default_y_key = 'z'
        self.default_z_key = 'rho'
        self.coords = ('r', 'z')

    def load_data(self):

        super().load_data()

        assert self.dim == 2, "Attempting to load 2D cylindrical data that is not 2D"

    @derived_field
    def cell_vol(self):

        dr, dz = self['dr']
        ones = np.ones_like(self['rho'])
        return (np.pi * (2*self['r'] - dr) * dr * dz * ones.T).T

    @derived_field
    def rtot(self):
        """ Array containing each cell's total radial distance from the origin. """

        return np.sqrt(self['r']**2 + self['z']**2)

    @derived_field
    def theta(self):
        """ Returns the polar angle corresponding to each cell. """

        return self['z'] / self.rtot

class Sedona3DCartPltFile(SedonaPltFileBase):

    def __init__(self, name, autoload=True, data=None):

        super().__init__(name, autoload, data)

        self.default_x_key = 'x'
        self.default_y_key = 'y'
        self.default_z_key = 'z'
        self.coords = ('x', 'y', 'z')

    def load_data(self):

        super().load_data()

        assert self.dim == 3, "Attempting to load 3D cartesian data that is not 3D"

    @derived_field
    def cell_vol(self):

        dx, dy, dz = self['dr']
        return dx * dy * dz * np.ones_like(self['rho'])

    @derived_field
    def rxy(self):
        """ Array containing each cell's total distance from the origin in the xy-plane. """

        return np.sqrt(self['x']**2 + self['y']**2)

    @derived_field
    def rtot(self):
        """ Array containing each cell's total radial distance from the origin. """

        return np.sqrt(self['x']**2 + self['y']**2 + self['z']**2)

    @derived_field
    def theta(self):
        """ Array containing the polar angle corresponding to each cell. """

        return self['z'] / self.rtot

    @derived_field
    def phi(self):
        """ Array containing the azimuthal angle corresponding the each cell. """

        return self['x'] / self.rxy