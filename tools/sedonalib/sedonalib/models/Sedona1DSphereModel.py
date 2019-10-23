from SedonaModel import SedonaBaseModel
import physical_constants as pc
import numpy as np

class Sedona1DSphereModel(SedonaBaseModel):


    def __init__(self,input_filename=None,dims=None,rout=None):

        # These variables will need to be defined
        # to constitute a legitimate model
        self.dims       = None
        self.n_zones    = 0
        self.dens       = None
        self.temp       = None
        self.erad       = None
        self.comp       = None
        self.r_edge     = None
        self.v_edge     = None
        self.elem_A     = []
        self.elem_Z     = []
        self.elem_list  = []
        self.time       = 1
        self.r_min      = 0
        self.type       = "standard"

        # read model from file
        if (input_filename is not None):
            fileformat = self.get_filename_format(input_filename)
            if (fileformat == "ascii"):
                self.read_ascii_file(input_filename)
            self.calculate_volumes()
        ######################

        # else return empty model
        if (dims is not None):
            self.n_zones = dims
            self.dims    = (dims,)
            self.dens    = np.zeros(self.n_zones)
            self.temp    = np.zeros(self.n_zones)
            self.erad    = np.zeros(self.n_zones)
            self.comp    = np.zeros((self.n_zones,0))
            self.v_edge  = np.zeros(self.n_zones)
            self.r_edge  = np.zeros(self.n_zones)
            self.vol     = np.zeros(self.n_zones)

            if (rout is not None):
                self.set_uniform_radial_edges(rout)

    def get_filename_format(self,name):
        ascii_extensions = [".dat",".mod",".txt"]
        hdf5_extensions = [".h5",".hdf5"]

        fileformat = ""
        if any(name.endswith(x) for x in ascii_extensions):
            fileformat = "ascii"
        if any(name.endswith(x) for x in hdf5_extensions):
            fileformat = "hdf5"
        return fileformat

    @property
    def r(self):
        r = np.zeros(self.n_zones)
        for i in range(self.n_zones):
            if (i==0):
                r[i] = 0.5*(self.r_min + self.r_edge[0])
            else:
                r[i] = 0.5*(self.r_edge[i] + self.r_edge[i-1])
        return r

    @property
    def v(self):
        v = np.zeros(self.n_zones)
        for i in range(self.n_zones):
            if (i==0):
                v[i] = 0.5*(0+ self.v_edge[0])
            else:
                v[i] = 0.5*(self.v_edge[i] + self.v_edge[i-1])
        return v


    @property
    def mass(self):
        return np.sum(self.vol*self.dens)

    @property
    def kinetic_energy(self):
        return np.sum(0.5*self.vol*self.dens*self.v_edge**2.0)

    def set_constant_velocity(self,v):
        self.v_edge.fill(v)


    def set_velocity(self,v):

        if (np.isscalar(v)):
            self.set_constant_velocity(v)
            return

        if (len(v) != self.n_zones):
            message = "length of array not equal to "
            message += "number of zones = " + str(self.n_zones) + "\n"
            raise ValueError(message)

        for i in range(self.n_zones):
            self.v_edge[i] = v[i]



    def set_uniform_radial_edges(self,rmax,rmin=0):

        dr = rmax/(1.0*self.n_zones)
        self.r_min = rmin
        self.r_edge = (np.arange(self.n_zones) + 1)*dr



    def calculate_volumes(self):

        self.vol = np.zeros(self.n_zones)
        for i in range(self.n_zones):
            if (i == 0):
                r0 = self.r_min
            else:
                r0 = self.r_edge[i-1]
            self.vol[i] = 4.0*np.pi/3.0*(self.r_edge[i]**3.0 - r0**3.0)


    def resize(self,new_nz):
        print "this is not implemented yet"

##############################################################
# Homologous profiles
##############################################################
    def set_homologous_constant_profile(self,mass,KE=None,time=86400.,vmax=None):

        self.time = time

        if (KE is None and vmax is None):
            raise exception("Must define either vmax or KE")

        if (vmax is None):
            vmax = (10.0/3.0*KE/mass)**(0.5)

        rho_t = mass/(4.0*pc.pi/3.0*vmax**3.0*self.time**3.0)

        self.r_min = 0.0
        dv = vmax/(1.0*self.n_zones)
        for i in range(self.n_zones):
            self.v_edge[i] = dv*(i+1.0)
            self.r_edge[i] = self.v_edge[i]*self.time
            self.dens[i] = rho_t

    def set_homologous_powerlaw_profile(self,mass,KE,t,n_out=10,n_in=1,vmax=None,min_density=1e-3):

        self.time = t
        v_t   = (2.0*KE/mass*(5.0 - n_in)*(n_out - 5.0)/(3.0 - n_in)/(n_out - 3.0))**0.5
        rho_t = (n_out - 3.0)*(3.0 - n_in)/(4.0*pc.pi)/(n_out - n_in)
        rho_t = rho_t*mass/(v_t**3*self.time**3.0)
        if (vmax is None):
            vmax = v_t*(min_density)**(-1.0/n_out)

        self.r_min = 0.0
        dv = vmax/(1.0*self.n_zones)
        for i in range(self.n_zones):
                self.v_edge[i] = dv*(i+1.0)
                self.r_edge[i] = self.v_edge[i]*self.time
                vm = dv*(i + 0.5)
                if (vm < v_t):
                    self.dens[i] = rho_t*(vm/v_t)**(-1.0*n_in)
                else:
                    self.dens[i] = rho_t*(vm/v_t)**(-1.0*n_out)


##############################################################
# Default supernova profiles
##############################################################
    def set_supernova_model(self,type,time=None,mass=None,KE=None,profile="powerlaw"):

        typeII_names = ["II","IIP","TYPEII","TYPEIIP","SNII","SNIIP"]
        typeIa_names = ["Ia","TYPEIa","SNIa"]

        if (type.upper() in typeII_names):
            if (mass is None):
                mass = 10.0*pc.m_sun
            if (KE is None):
                KE = 1e51

        elif (type.upper() in typeIa_names):
            if (mass is None):
                mass = 1.0*pc.m_sun


##############################################################
# Reading in functions
##############################################################

    def read_ascii_file(self,name):

        with open(name) as fin:

            # read header
            geometry,rtype = fin.readline().split()
            nz, rin, t, nelem  = fin.readline().split()
            self.n_zones = int(nz)
            self.time    = float(t)
            n_elements   = int(nelem)
            self.r_min   = float(rin)

            # determine model type
            if (rtype == "SNR" or rtype == "homologous"):
                self.r_min = self.r_min*self.time
                self.type = "homologous"
            elif (rtype == "standard"):
                self.type = "standard"
            else:
                message = "Unknown model type " + rtype
                message += " in header\n"
                raise ValueError(message)

            # read in element list
            self.elem_list = fin.readline().split()
            if (len(self.elem_list) != n_elements):
                message = "Badly formatted sedona model, "
                message += "n_elements in header not equal to list"
                raise ValueError(message)

            # read the rest of the data
            data = np.genfromtxt(fin)

        # close the file

        if self.n_zones != len(data):
            raise ValueError('badly formatted sedona model: nzones != number of data rows')

        # determine elements
        self.elem_Z = np.zeros(n_elements,dtype='i')
        self.elem_A = np.zeros(n_elements,dtype='i')
        for i,e in enumerate(self.elem_list):
            Z,A =  e.split(".")
            self.elem_Z[i] = int(Z)
            self.elem_A[i] = int(A)

        if (self.type == "standard"):
            self.r_edge = data[:, 0]
            self.v_edge = data[:, 1]
            self.dens   = data[:, 2]
            self.temp   = data[:, 3]
            self.comp   = data[:,4:]
        elif (self.type == "SNR" or type == "homologous"):
            self.v_edge = data[:, 0]
            self.dens   = data[:, 1]
            self.temp   = data[:, 2]
            self.comp   = data[:,3:]
            self.r_edge = self.vel*self.time
        else:
            raise Exception("Unknown model type " + type)



##############################################################
# Writing out functions
##############################################################

    def check_model_validity(self):

        if (self.check_base_model_validity() is False):
            return False

        if (self.n_zones <= 0):
            return False

        if (len(self.elem_Z) <= 0):
            print 'Error: no elements defined in composition'
            return False

        if (self.dens is None):
            print 'Error: density has not been set'
            return False

        return True


    def normalize_composition(self):

        nel = len(self.elem_A)
        for i in range(self.n_zones):
            norm = np.sum(self.comp[i,:])

            if (norm == 0):
                for j in range(nel):
                    self.comp[i,j] = 1.0/(1.0*nel)
            else:
                self.comp[i,:] /= norm



    def write(self,outfile):

        valid = self.check_model_validity()
        if (not valid):
            message = "Model is not valid or complete\n"
            raise Exception(message)

        self.normalize_composition()

        fileformat = self.get_filename_format(outfile)
        if (fileformat == "ascii"):
            self.write_ascii_file(outfile)
        if (fileformat == "hdf5"):
            self.write_hdf5_file(outfile)

        # through error if file format unknown
        if (fileformat == ""):
            message = "Unknown filename extension, allowed values are "
            message += str(ascii_extensions + hdf5_extensions)
            raise Exception(message)


    def write_hdf5_file(self,outfile):

        import h5py
        fout = h5py.File(outfile, 'w')
        fout.create_dataset('r_out',data=self.r_edge,dtype='d')
        fout.create_dataset('r_min',data=self.r_min,dtype='d')

        fout.create_dataset('time',data=[self.time],dtype='d')
        fout.create_dataset('Z',data=self.elem_Z,dtype='i')
        fout.create_dataset('A',data=self.elem_A,dtype='i')
        fout.create_dataset('rho',data=self.dens,dtype='d')
        fout.create_dataset('temp',data=self.temp,dtype='d')
        fout.create_dataset('v',data=self.v_edge,dtype='d')
        fout.create_dataset('erad',data=self.erad,dtype='d')
        fout.create_dataset('comp',data=self.comp,dtype='d')
        fout.close()


    def write_ascii_file(self,outfile):

        with open(outfile, 'w') as f:

            # write header
            n_elems = len(self.elem_Z)
            f.write("1D_sphere " + self.type + "\n")
            f.write("%d %f %f %d\n" % (self.n_zones,
                                       self.r_min,
                                       self.time,
                                       n_elems))

            # write element list
            for i in range(n_elems):
                f.write("%d.%d " %(self.elem_Z[i],self.elem_A[i]))
            f.write("\n")

            # write zone data
            for i in range(self.n_zones):

                line = ""
                if (self.type != "homologous"):
                    line = ['%e' % ff for ff in [self.r_edge[i]]]
                line += ['%e' % ff for ff in [self.v_edge[i], self.dens[i], self.temp[i]]]
                line += ['%e' % ff for ff in self.comp[i].tolist()]
                f.write(' '.join(line) + '\n')



##############################################################
# Print out and plot out functions
##############################################################

    def __str__(self):

        lbreak = "---------------\n"

        ret_str = lbreak
        ret_str += "Sedona model\n" + lbreak
        ret_str += "geometry   = 1D_sphere\n"
        ret_str += "type       = " + self.type + "\n"
        ret_str += "n_zones    = {:0}\n".format(self.n_zones)
        ret_str += "n_elements = {:0}\n".format(len(self.elem_Z))
        ret_str += lbreak

        for i in range(len(self.elem_Z)):
            ret_str += " {:0}.{:1}\n".format(self.elem_Z[i],self.elem_A[i])

        ret_str += lbreak
        ret_str += "mass = {:0}\n".format(self.mass)

        return ret_str



    def plot(self, show=True):

        try:
            import matplotlib.pyplot as plt
            from matplotlib.pyplot import cm
        except Exception as e:
            raise e

        try:
            import seaborn as sns
        except:
            sb = False
        else:
            sb = True

        if sb:
            sns.set_style('ticks')

        fig, axarr = plt.subplots(nrows=2,ncols=2,figsize=(12,9))

        colors = ['b','g','r','purple','c','m','y','k','orange','indigo','violet']
        ls = ['-', '--', ':', '-.']
        styles = [{'color':c, 'ls':linestyle} for c in colors for linestyle in ls]

        x = self.v_edge
        xlabel = 'velocity (cm/s)'

        # plot density structure
        axarr[0,0].semilogy(x,self.dens)
        axarr[0,0].set_ylabel('rho (g / cm3)')
        axarr[0,0].set_xlabel(xlabel)

        # plot temperature structure
        axarr[0,1].plot(x, self.temp)
        axarr[0,1].set_ylabel('T (Kelvin)')
        axarr[0,1].set_xlabel(xlabel)

        # plot velocity
        axarr[1,0].plot(self.r_edge, self.v_edge)
        axarr[1,0].set_ylabel('velocity (cm/s)')
        axarr[1,0].set_xlabel("radius (cm)")

        # plot mass fractions structure
        for i,e in enumerate(self.elem_list):
            axarr[1,1].semilogy(x,self.comp[:,i])
        axarr[1,1].set_ylabel('mass fraction')
        axarr[1,1].set_ylim(1e-6,1.5)
        axarr[1,1].set_xlabel(xlabel)


        for ax in axarr.ravel():
            ax.minorticks_on()

        if sb:
            sns.despine()

#        fig.tight_layout()
#        fig.subplots_adjust(top=0.8)
        handles, labels = axarr[1,1].get_legend_handles_labels()
        fig.legend(handles,labels, loc='lower right', ncol=5)

        if show:
            plt.ion()
            fig.show()
            j = raw_input("press any key to continue>")

        return (fig, axarr)
