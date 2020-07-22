from .SedonaModel import SedonaBaseModel
#from Sedona1DSphereModel import Sedona1DSphereModel
import sedonalib.physical_constants as pc
import numpy as np
import h5py

class Sedona2DCylnModel(SedonaBaseModel):


    def __init__(self,input_filename=None,dims=None,rout=None):

        # These variables will need to be defined
        # to constitute a legitimate model
        self.dims       = None
        self.dens       = None
        self.temp       = None
        self.erad       = None
        self.comp       = None
        self.vx         = None
        self.vz         = None
        self.elem_A     = []
        self.elem_Z     = []
        self.elem_list  = []
        self.time       = 1
        self.dr         = None

        if (dims is not None):
            self.set_dims(dims)

        # read model from file
        if (input_filename is not None):
            self.read(input_filename)


    ###############################################
    # Read model file
    ###############################################

    def read(self,infile):


        fin = h5py.File(infile, 'r')
        self.vx = np.array(fin['vx'],dtype='d')
        self.vz = np.array(fin['vz'],dtype='d')
        self.dens = np.array(fin['rho'],dtype='d')
        self.temp = np.array(fin['temp'],dtype='d')
        self.comp = np.array(fin['comp'],dtype='d')
        self.erad = np.array(fin['erad'],dtype='d')
        self.elem_A = np.array(fin['A'],dtype='i')
        self.elem_Z = np.array(fin['Z'],dtype='i')
        self.dr = np.array(fin['dr'],dtype='d')
        self.time = float((fin['time'])[0])
        self.dims = self.dens.shape


    ###############################################
    # Return basic properties
    ###############################################
    @property
    def mass(self):

        x = self.x
        mass = 0.0

        for i in range(self.dims[0]):
            for j in range(self.dims[1]):
                vol  = 2.0*np.pi*x[i]*self.dr[0]*self.dr[1]
                mass += vol*self.dens[i,j]
        return mass

    @property
    def kinetic_energy(self):

        x = self.x
        ke = 0.0

        for i in range(self.dims[0]):
            for j in range(self.dims[1]):
                vol  = 2.0*np.pi*x[i]*self.dr[0]*self.dr[1]
                vsq = self.vx[i,j]**2.0 + self.vz[i,j]**2.0
                ke += 0.5*vol*self.dens[i,j]*vsq

        return ke

    @property
    def vtot(self):
        return (self.vx**2.0 + self.vz**2.0)**0.5

    @property
    def x(self):
        return (np.arange(self.dims[0])+0.5)*self.dr[0]

    @property
    def z(self):
        zmin = -1.0*self.dims[1]*self.dr[1]/2.0
        return zmin + (np.arange(self.dims[1]) + 0.5)*self.dr[1]

    @property
    def theta(self):
        z = self.z
        x = self.x
        theta = np.zeros(self.dims)
        for i in range(len(x)):
            for j in range(len(z)):
                r = (x[i]**2.0 + z[j]**2.0)**0.5
                theta[i,j] = np.arccos(z[j]/r)
        return theta

    @property
    def r(self):
        z = self.z
        x = self.x
        r = np.zeros(self.dims)
        for i in range(len(x)):
            for j in range(len(z)):
                r[i,j] = (x[i]**2.0 + z[j]**2.0)**0.5
        return r

    def get_spherical_coordinates(self):
        z = self.z
        x = self.x
        r = np.zeros(self.dims)
        theta = np.zeros(self.dims)
        for i in range(len(x)):
            for j in range(len(z)):
                r[i,j] = (x[i]**2.0 + z[j]**2.0)**0.5
                theta[i,j] = np.arccos(z[j]/r[i,j])
        return r,theta

    ###############################################
    # Set basic print_properties
    ###############################################

    def set_grid(self,dims,dr):

        self.set_dims(dims)
        self.set_dr(dr)

    def set_dims(self, dims):

        if (np.isscalar(dims)):
            self.dims = (dims,2*dims)
        else:
            self.dims = dims

        self.vx = np.zeros(self.dims)
        self.vz = np.zeros(self.dims)
        self.dens = np.zeros(self.dims)
        self.temp = np.zeros(self.dims)
        self.erad = np.zeros(self.dims)

        if (self.elem_A is None):
            nelem = 0
        else:
            nelem = len(self.elem_A)
        self.comp = np.zeros(self.dims +(nelem,))

    def set_dr(self, d):

        self.dr = np.zeros(2)
        if (np.isscalar(d)):
            self.dr[0] = d
            self.dr[1] = d
        else:
            self.dr[0] = d[0]
            self.dr[1] = d[1]

    ###############################################
    # Remap
    ###############################################

    def set_homologous_profile(self,mass,KE,type="powerlaw",time=86400,vmax=None,n_out=10,n_in=1,min_density=1e-3):

        from Sedona1DSphereModel import Sedona1DSphereModel

        mod1d = Sedona1DSphereModel(dims=self.dims[0])
        mod1d.set_homologous_profile(mass,KE,type=type,vmax=vmax,n_out=n_out,n_in=n_in,min_density=min_density)
        self.remap_from(mod1d)


    def remap_from(self,mod):
        self.remap_from_1DSphereModel(mod)

    def remap_from_1DSphereModel(self,mod1d):

        rho_fac = 1e-10

        self.time = mod1d.time

        # copy over elements
        self.elem_A = list(mod1d.elem_A)
        self.elem_Z = list(mod1d.elem_Z)
        nelem = len(self.elem_A)

        # reset dimensions if not set already
        n1d = (mod1d.dims)[0]
        rmax_1d = mod1d.rmax
        rmin_1d = mod1d.rmin
        if (self.dims is None):
            self.set_dims(mod1d.dims)
            dr = rmax_1d/(1.0*n1d)
            self.set_dr(dr)

        if (self.dr is None):

            dr = rmax_1d/(1.0*n1d)
            self.set_dr(dr)

        r1d = mod1d.r
        x = self.x
        z = self.z

        for i in range(self.dims[0]):
            for j in range(self.dims[1]):

                r = (x[i]**2 + z[j]**2.0)**0.5
                if (r < rmax_1d and r > rmin_1d):
                    self.temp[i,j] = (np.interp([r],r1d,mod1d.temp))[0]
                    self.dens[i,j] = (np.interp([r],r1d,mod1d.dens))[0]

                    if (mod1d.erad is None):
                        self.erad[i,j] = 0.0
                    else:
                        self.erad[i,j] = (np.interp([r],r1d,mod1d.erad))[0]

                    v = (np.interp([r],r1d,mod1d.v))[0]
                    self.vx[i,j]  = v*x[i]/r
                    self.vz[i,j]  = v*z[j]/r
                    for l in range(nelem):
                        self.comp[i,j,l] = (np.interp([r],r1d,mod1d.comp[:,l]))[0]

                elif (r >= rmax_1d):
                    self.temp[i,j] = mod1d.temp[-1]
                    self.dens[i,j] = mod1d.dens[-1]*rho_fac
                    self.vx[i,j]   = 0
                    self.vz[i,j]   = 0.0
                    self.erad[i,j] = 0.0
                    for l in range(nelem):
                        self.comp[i,j,l] = mod1d.comp[-1,l]

                else:
                    self.temp[i,j] = mod1d.temp[0]
                    self.dens[i,j] = mod1d.dens[0]*rho_fac
                    self.vx[i,j]   = 0
                    self.vz[i,j]   = 0.0
                    self.erad[i,j] = 0.0
                    for l in range(nelem):
                        self.comp[i,j,l] = mod1d.comp[0,l]


    ###############################################
    # Plot, Print and Write
    ###############################################

    def __str__(self):

        lbreak = "---------------\n"

        ret_str = lbreak
        ret_str += "Sedona model\n" + lbreak
        ret_str += "geometry   = 2D_cyln\n"
        ret_str += "dims       = ({:0},{:1})\n".format(self.dims[0],self.dims[1])
        nzones = self.dims[0]*self.dims[1]
        ret_str += "n_zones    = {:0}\n".format(nzones)
        ret_str += "n_elements = {:0}\n".format(len(self.elem_Z))
        ret_str += "dr         = ({:0},{:1})\n".format(self.dr[0],self.dr[1])
        ret_str += lbreak

        for i in range(len(self.elem_Z)):
            ret_str += " {:0}.{:1}\n".format(self.elem_Z[i],self.elem_A[i])

        ret_str += lbreak
        ret_str += "mass = {:0} g\n".format(self.mass)
        ret_str += "kinetic energy = {:0} erg\n".format(self.kinetic_energy)

        return ret_str


    def check_model_validity(self):

        if (self.check_base_model_validity() is False):
            return False

        if (self.dr is None):
            print ('Error: size of zones dr not set')
            return False

        if (self.dims is None):
            print ('Error: dimensions of grid not set')
            return False

        return True

    def write(self, outfile):

        if (not self.check_model_validity()):
            message = "Model is not valid or complete\n"
            raise Exception(message)

        fout = h5py.File(outfile, 'w')
        fout.create_dataset('dr',data=self.dr,dtype='d')
        fout.create_dataset('time',data=[self.time],dtype='d')
        fout.create_dataset('Z',data=self.elem_Z,dtype='i')
        fout.create_dataset('A',data=self.elem_A,dtype='i')
        fout.create_dataset('rho',data=self.dens,dtype='d')
        fout.create_dataset('temp',data=self.temp,dtype='d')
        fout.create_dataset('vx',data=self.vx,dtype='d')
        fout.create_dataset('vz',data=self.vz,dtype='d')
        fout.create_dataset('erad',data=self.erad,dtype='d')
        fout.create_dataset('comp',data=self.comp,dtype='d')
        fout.close()


    def plot(self,log=True):

        try:
            import matplotlib.pyplot as plt
            from matplotlib.pyplot import cm
        except Exception as e:
            raise e

        fig,axs = plt.subplots(1,2)
        plt.ion()
        im1 = axs[0].matshow(np.swapaxes(np.log10(self.dens),0,1))
        axs[0].set_title("density")
        fig.colorbar(im1, ax=axs[0])
        im2 = axs[1].matshow(np.swapaxes(self.temp,0,1))
        axs[1].set_title("temperature")
        fig.colorbar(im2, ax=axs[1])
        plt.show()
        self.wait_for_key()

        plt.ion()
        fig,axs = plt.subplots(1,2)
        im1 = axs[0].matshow(np.swapaxes(self.vx,0,1))
        axs[0].set_title("x-velocity")
        fig.colorbar(im1, ax=axs[0])
        im2 = axs[1].matshow(np.swapaxes(self.vz,0,1))
        axs[1].set_title("z-velocity")
        fig.colorbar(im2, ax=axs[1])
        plt.show()
        self.wait_for_key()



        nelem = len(self.elem_Z)
        for i in range(0,nelem,2):

            plt.clf()
            fig,axs = plt.subplots(1,2)

            im1 = axs[0].matshow(np.swapaxes(self.comp[:,:,i],0,1))
            axs[0].set_title("mass fraction Z = " + str(self.elem_Z[i]))
            fig.colorbar(im1, ax=axs[0])
            if (i+1 < nelem):
                arr = np.swapaxes(self.comp[:,:,i+1],0,1)
                if (log):
                    arr = np.log10(arr)
                im2 = axs[1].matshow(arr)
                axs[1].set_title("mass fraction Z = " + str(self.elem_Z[i+1]))
                fig.colorbar(im2, ax=axs[1])
            plt.show()
            self.wait_for_key()
