from SedonaModel import SedonaBaseModel
import physical_constants as pc
import numpy as np
import h5py

class Sedona3DCartModel(SedonaBaseModel):


    def __init__(self,input_filename=None,dims=None,rout=None):

        # These variables will need to be defined
        # to constitute a legitimate model
        self.dims       = dims
        self.dens       = None
        self.temp       = None
        self.erad       = None
        self.comp       = None
        self.vx         = None
        self.vy         = None
        self.vz         = None
        self.elem_A     = []
        self.elem_Z     = []
        self.elem_list  = []
        self.time       = 1
        self.rmin       = None
        self.dr         = None

        if (dims is not None):
            self.set_dims(dims)

        # read model from file
        if (input_filename is not None):
            self.read(input_filename)
        ######################

    ###############################################
    # Read model file
    ###############################################

    def read(self, infile):

        fin = h5py.File(infile, 'r')
        self.vx = np.array(fin['vx'],dtype='d')
        self.vy = np.array(fin['vy'],dtype='d')
        self.vz = np.array(fin['vz'],dtype='d')
        self.dens = np.array(fin['rho'],dtype='d')
        self.temp = np.array(fin['temp'],dtype='d')
        self.comp = np.array(fin['comp'],dtype='d')
        self.erad = np.array(fin['erad'],dtype='d')
        self.elem_A = np.array(fin['A'],dtype='i')
        self.elem_Z = np.array(fin['Z'],dtype='i')
        self.dr = np.array(fin['dr'],dtype='d')
        self.rmin = np.array(fin['rmin'],dtype='d')
        self.time = float((fin['time'])[0])
        self.dims = self.dens.shape


    ###############################################
    # Return basic properties
    ###############################################
    @property
    def mass(self):
        vol  = self.dr[0]*self.dr[1]*self.dr[2]
        return np.sum(vol*self.dens)

    @property
    def kinetic_energy(self):
        vol  = self.dr[0]*self.dr[1]*self.dr[2]
        vtotsq = (self.vx**2.0 + self.vy**2.0 + self.vz**2.0)
        return np.sum(0.5*vol*self.dens*vtotsq)

    @property
    def vtot(self):
        return (self.vx**2.0 + self.vy**2.0 + self.vz**2.0)**0.5

    @property
    def x(self):
        return self.rmin[0] + np.arange(self.dims[0])*self.dr[0]
    @property
    def y(self):
        return self.rmin[1] + np.arange(self.dims[1])*self.dr[1]
    @property
    def z(self):
        return self.rmin[2] + np.arange(self.dims[2])*self.dr[2]

    ###############################################
    # Set basic properties
    ###############################################

    def set_grid(self,dims,dr,rmin):

        self.set_dims(dims)
        self.set_dr(dr)
        self.set_rmin(rmin)


    def set_dims(self, dims):

        if (np.isscalar(dims)):
            self.dims = (dims,dims,dims)
        else:
            self.dims = dims

        self.vx = np.zeros(self.dims)
        self.vy = np.zeros(self.dims)
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

        self.dr = np.zeros(3)
        if (np.isscalar(d)):
            self.dr[0] = d
            self.dr[1] = d
            self.dr[2] = d
        else:
            self.dr[0] = d[0]
            self.dr[1] = d[1]
            self.dr[2] = d[2]


    def set_rmin(self,r):

        self.rmin = np.zeros(3)
        if (np.isscalar(r)):
            self.rmin[0] = r
            self.rmin[1] = r
            self.rmin[2] = r
        else:
            self.rmin[0] = r[0]
            self.rmin[1] = r[1]
            self.rmin[2] = r[2]

    def set_vx(self,v):

        if (np.isscalar(v)):
            self.vx.fill(v)
            return

        if (v.shape != self.dims):
            message =  "dimensions of array not equal to "
            message += "grid dimensions = " + str(self.dims) + "\n"
            raise ValueError(message)

        self.vx = np.copy(v)

    def set_vy(self,v):

        if (np.isscalar(v)):
            self.vy.fill(v)
            return

        if (v.shape != self.dims):
            message =  "dimensions of array not equal to "
            message += "grid dimensions = " + str(self.dims) + "\n"
            raise ValueError(message)

        self.vy = np.copy(v)

    def set_vz(self,v):

        if (np.isscalar(v)):
            self.vz.fill(v)
            return

        if (v.shape != self.dims):
            message =  "dimensions of array not equal to "
            message += "grid dimensions = " + str(self.dims) + "\n"
            raise ValueError(message)

        self.vz = np.copy(v)

    def set_velocity(self,vx,vy,vz):
        self.set_vx(vx)
        self.set_vy(vy)
        self.set_vz(vz)

    def set_homologous_velocities(self,time=None):

        if (time is not None):
            self.time = time

        for i in range(self.dims[0]):
            for j in range(self.dims[1]):
                for k in range(self.dims[2]):
                    x = self.rmin[0] + self.dr[0]*(i + 0.5)
                    y = self.rmin[1] + self.dr[1]*(j + 0.5)
                    z = self.rmin[2] + self.dr[2]*(k + 0.5)
                    self.vx[i,j,k] = x/self.time
                    self.vy[i,j,k] = y/self.time
                    self.vz[i,j,k] = z/self.time


    def check_homology(self):

        try:
            import matplotlib.pyplot as plt
        except Exception as e:
            raise e

        plt.plot(self.x,self.vx[:,self.dims[1]/2,self.dims[2]/2],'o')
        plt.plot(self.y,self.vy[self.dims[0]/2,:,self.dims[2]/2],'o')
        plt.plot(self.z,self.vz[self.dims[0]/2,self.dims[1]/2,:],'o')
        plt.plot(self.x,self.x/self.time)

        plt.legend(['X','Y','Z'])
        plt.ion()
        plt.show()
        self.wait_for_key()

    def homologously_evolve(self,tnew,told=None):

        if (told is not None):
            self.time = told

        fac = tnew/self.time

        for i in range(3):
            self.dr[i]   = self.dr[i]*fac
            self.rmin[i] = self.rmin[i]*fac

        self.dens *= 1.0/fac**3.0
        self.temp *= 1.0/fac
        self.time = tnew


    ###############################################
    # Plot, Print and Write
    ###############################################

    def check_model_validity(self):

        if (self.check_base_model_validity() is False):
            return False

        if (self.dr is None):
            print 'Error: size of zones dr not set'
            return False

        if (self.dims is None):
            print 'Error: dimensions of grid not set'
            return False

        if (self.rmin is None):
            print 'Error: rmin boundary of grid not set'
            return False

        return True

    def write(self, outfile):

        if (not self.check_model_validity()):
            message = "Model is not valid or complete\n"
            raise Exception(message)

        fout = h5py.File(outfile, 'w')
        fout.create_dataset('dr',data=self.dr,dtype='d')
        fout.create_dataset('rmin',data=self.rmin,dtype='d')

        fout.create_dataset('time',data=[self.time],dtype='d')
        fout.create_dataset('Z',data=self.elem_Z,dtype='i')
        fout.create_dataset('A',data=self.elem_A,dtype='i')
        fout.create_dataset('rho',data=self.dens,dtype='d')
        fout.create_dataset('temp',data=self.temp,dtype='d')
        fout.create_dataset('vx',data=self.vx,dtype='d')
        fout.create_dataset('vy',data=self.vy,dtype='d')
        fout.create_dataset('vz',data=self.vz,dtype='d')
        fout.create_dataset('erad',data=self.erad,dtype='d')
        fout.create_dataset('comp',data=self.comp,dtype='d')
        fout.close()



    def plot(self):

        try:
            import matplotlib.pyplot as plt
            from matplotlib.pyplot import cm
#            from mayavi import mlab
        except Exception as e:
            raise e

        nxm = int(self.dims[0]/2)
        nym = int(self.dims[1]/2)
        nzm = int(self.dims[2]/2)

        fix,axs = plt.subplots(2,3)
        axs[0,0].matshow(np.log10(self.density[:,:,nzm]))
    #    axs[0,0].colorbar()
        axs[0,1].matshow(self.temp[:,:,nzm])
        axs[0,2].matshow(self.vtot[:,:,nzm])
        axs[1,0].matshow(self.vx[:,:,nzm])
        axs[1,1].matshow(self.vy[:,:,nzm])
        axs[1,2].matshow(self.vz[:,:,nzm])

#        plt.colorbar()

#        mlab.figure(1, fgcolor=(1, 1, 1), bgcolor=(1, 1, 1))
#        mlab.contour3d(np.log10(self.dens + 1e-20),contours=20,opacity=0.2,colormap='jet' )

        plt.ion()
        plt.show()
        self.wait_for_key()
