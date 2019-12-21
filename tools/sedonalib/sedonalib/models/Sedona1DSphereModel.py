from SedonaModel import SedonaBaseModel
from Sedona2DCylnModel import Sedona2DCylnModel
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
        self.m_edge     = None
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
            self.set_dims(dims)

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
    def rmax(self):
        return self.r_edge[-1]
    @property
    def rmin(self):
        return self.r_min

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
        self.calculate_volumes()
        return np.sum(self.vol*self.dens)

    @property
    def kinetic_energy(self):
        return np.sum(0.5*self.vol*self.dens*self.v_edge**2.0)

    def set_constant_velocity(self,v):
        self.v_edge.fill(v)


    def set_dims(self,dims):
        self.n_zones = dims
        self.dims    = (dims,)
        self.dens    = np.zeros(self.n_zones)
        self.temp    = np.zeros(self.n_zones)
        self.erad    = np.zeros(self.n_zones)
        self.comp    = np.zeros((self.n_zones,0))
        self.v_edge  = np.zeros(self.n_zones)
        self.r_edge  = np.zeros(self.n_zones)
        self.m_edge  = np.zeros(self.n_zones)
        self.vol     = np.zeros(self.n_zones)


    def set_grid(self,dims,rmax=None,rmin = 0.0):
        self.set_dims(dims)
        if (rmax is not None):
            self.set_uniform_radial_edges(rmax,rmin=rmin)

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

        self.vol    = np.zeros(self.n_zones)
        self.m_edge = np.zeros(self.n_zones)
        for i in range(self.n_zones):
            if (i == 0):
                r0 = self.r_min
                m0 = 0.0
            else:
                r0 = self.r_edge[i-1]
                m0 = self.m_edge[i-1]
            self.vol[i] = 4.0*np.pi/3.0*(self.r_edge[i]**3.0 - r0**3.0)
            self.m_edge[i] = m0 + self.vol[i]*self.rho[i]

    def resize(self,new_nz):
        print "this is not implemented yet"

##############################################################
# Set composition functions
##############################################################
    def set_shell_composition(self,comp,mrange=None,vrange=None,rrange=None):

        err = False
        if (mrange is not None and vrange is not None): err = True
        if (mrange is not None and rrange is not None): err = True
        if (vrange is not None and rrange is not None): err = True
        if (err):
            mess = "Can define only one of mrange, vrange, or rrange"
            raise ValueError(mess)

        # set composition in mass range
        if (mrange is not None):

            if (len(mrange) != 2):
                raise ValueError("mrange must have format [mass_in,mass_out]")
            if (mrange[0] > mrange[1]):
                raise ValueError("mrange must have format [mass_in,mass_out]")

            # make sure mass coordinates are calculated
            self.calculate_volumes()

            b = (self.m_edge >= mrange[0])*(self.m_edge <= mrange[1])
            self.comp[b,:] = comp

        # set composition in radial range
        if (rrange is not None):

            if (len(rrange) != 2):
                raise ValueError("rrange must have format [r_in,r_out]")
            if (rrange[0] > rrange[1]):
                raise ValueError("rrange must have format [r_in,r_out]")

            b = (self.r_edge >= rrange[0])*(self.r_edge <= rrange[1])
            self.comp[b,:] = comp

        # set composition in velocity range
        if (vrange is not None):

            if (len(vrange) != 2):
                raise ValueError("rrange must have format [v_in,v_out]")
            if (vrange[0] > vrange[1]):
                raise ValueError("rrange must have format [v_in,v_out]")
            # check if velocity is monotonic
            for i in range(self.v_edge.size-1):
                if (self.v_edge[i+1] < self.v_edge[i]) :
                    raise ValueError("Velocity must be monotonic to set shell")

            b = (self.v_edge >= vrange[0])*(self.v_edge <= vrange[1])
            self.comp[b,:] = comp


##############################################################
# Homologous density profiles
##############################################################

    def set_homologous_profile(self,mass,KE,type="powerlaw",time=86400,vmax=None,n_out=10,n_in=1,min_density=1e-3):

        if (self.dims is None):
            message = "Error == Must set model dimensions before setting homologous profile"
            raise ValueError(message)

        if (type == "constant"):
            self.set_homologous_constant_profile(mass,KE=KE,time=time,vmax=vmax)
        elif (type == "powerlaw"):
            self.set_homologous_powerlaw_profile(mass,KE,time=time,n_out=n_out,n_in=n_in,vmax=vmax,min_density=min_density)
        elif (type == "exponential"):
            self.set_homologous_exponential_profile(mass,KE,time=time,vmax=vmax)
        else:
            message = "Error == Unknown homolgous profile type " + type
            raise ValueError(message)



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

        self.calculate_volumes()


    def set_homologous_powerlaw_profile(self,mass,KE,time=86400.,n_out=10,n_in=1,vmax=None,min_density=1e-3):

        self.time = time
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

        self.calculate_volumes()


    def set_homologous_exponential_profile(self,mass,KE,time=86400.,vmax=None,min_density=1e-3):

        self.time = time
        v_e   = (KE/6.0/mass)**(0.5)
        rho0  = mass/8.0/np.pi/(v_e*self.time)**3.0

        if (vmax is None):
            vmax = v_e*np.log(1.0/min_density)

        self.r_min = 0.0
        dv = vmax/(1.0*self.n_zones)
        for i in range(self.n_zones):
            vm = dv*(i + 0.5)
            self.v_edge[i] = dv*(i+1.0)
            self.r_edge[i] = self.v_edge[i]*self.time
            self.dens[i] = rho0*np.exp(-vm/v_e)

        self.calculate_volumes()


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
# Remapping functions
##############################################################
    def remap_to(self,type):

        if (type == "2D_cyln"):
            newmod = Sedona2DCylnModel()
            newmod.remap_from(self)
            return newmod

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

        self.dims = self.n_zones
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

        self.calculate_volumes()
        for i in range(len(self.elem_Z)):
            elem_mass = np.sum(self.vol*self.dens*self.comp[:,i])
            ret_str += " {0:3.0f}.{1:}  ".format(self.elem_Z[i],self.elem_A[i])
            ret_str += " {0:.4e} ({1:.4e} msun)\n".format(elem_mass,elem_mass/pc.m_sun)
        ret_str += lbreak
        ret_str += "mass = {0:.4e} g ({1:.4e} msun)\n".format(self.mass,self.mass/pc.m_sun)
        ret_str += "kinetic energy = {0:.4e} erg\n".format(self.kinetic_energy)

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
