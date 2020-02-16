import numpy as np
import h5py

class SedonaPltFile():


    def __init__(self,name):

        self.name = name
        self.load_data()


    def __getitem__(self, key):
        return self.data[key]

    def load_data(self):

        ds         = h5py.File(self.name)
        keys       = ds.keys()
        self.data  = {key: np.array(ds[key], dtype = np.float64) for key in keys}
        ds.close()

        self.data_status 		= "loaded"
        self.data_keys 			= self.data.keys()
        self.dim 				= len(self.data['rho'].shape)

        print("Found the modelfile to be a ", self.dim, "dimensional setup.")

        if self.dim == 1:
            self.default_x_key = 'r'
            self.default_y_key = 'rho'
        elif self.dim == 2:
            self.default_x_key = 'r'
            self.default_y_key = 'z'
            self.default_z_key = 'rho'
        elif self.dim == 3:
            self.default_x_key = 'x'
            self.default_y_key = 'y'
            self.default_z_key = 'z'
        else:
            print("Not setting any default options for the plot")