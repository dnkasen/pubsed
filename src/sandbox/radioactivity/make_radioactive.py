import h5py 

# python script defining radioactive isotope data and writing to either a
# a hdf5 or a ascii file

#
class isotope:
	name         = ''      # name of isotope in the format Z.A
	daughter     = ''      # name of daughter istope in Z.A
	t_decay      = 0.0     # decay time in seconds
	E_decay      = 0.0     # average energy per decay in MeV
	p_frac       = 0.0     # percent of decay energy that goes to positron kinetic 
	g_energy     = []      # list of energies of emitted gamma rays in MeV
	g_branch     = []      # probability of each gamma-ray line

	def add_to_hdf5(self,fname):

		f = h5py.File(fname,'a')
		base = self.name + "/"
		grp = f.create_group(base)
		f.create_dataset(base + "t_decay",data = self.t_decay)
		f.create_dataset(base + "daughter",data = self.daughter)
		f.create_dataset(base + "energy",data=self.E_decay)
		f.create_dataset(base + "positron_fraction",data=self.p_frac)
		f.create_dataset(base + "gamma_branching",data=self.g_branch)
		f.create_dataset(base + "gamma_energies",data=self.g_energy)
		f.close()

	def add_to_ascii(self,fname):
		# need to add code in some format

##################################

# write out some iosotoeps
fname = "radioactive_data.h5"
f = h5py.File(fname,'w')
f.close()

iso = isotope()

### Ni56
iso.name      = "28.56"
iso.daughter  = "27.56"
iso.t_decay   = 757241.7
iso.E_decay   = 1.728
iso.p_frac    = 0.0
iso.g_energy  = [0.15838, 0.2695, 0.48044, 0.74995, 0.81185, 1.5618]
iso.g_branch  = [0.988, 0.365, 0.365, 0.495, 0.86, 0.14]
iso.add_to_hdf5(fname)

### Co56
iso.name      = "27.56"
iso.daughter  = "26.56"
iso.t_decay   = 9627378.
iso.E_decay   = 3.73
iso.p_frac    = 0.0321
iso.g_energy  = [0.39, 0.00191, 0.00311, 0.999399, 0.00049, 
                              0.00073, 0.01421, 0.00111, 0.1405, 0.00055, 
                              0.00132, 0.00094, 0.02252, 0.00049, 0.6646, 
                              0.001224, 0.04283, 0.0018, 0.00074, 0.000616, 
                              0.1541, 0.0064, 0.00707, 0.03016, 0.0777, 
                              0.00377, 0.00388, 0.00118, 0.0008, 0.00059, 
                              0.1697, 0.00019, 0.01036, 0.03209, 0.07923, 
                              0.018759, 0.00949, 0.001955, 0.000167]

iso.g_branch  = [0.511, 0.733514, 0.787743, 0.84677, 0.852732, 
                              0.89651, 0.977372, 0.996948, 1.037843, 1.088894, 
                              1.140368, 1.159944, 1.175101, 1.198888, 1.238288, 
                              1.3354, 1.360212, 1.442746, 1.462322, 1.640475, 
                              1.771357, 1.810757, 1.963741, 2.015215, 2.034791, 
                              2.113135, 2.212944, 2.276131, 2.37324, 2.52309, 
                              2.5985, 2.657527, 3.009645, 3.202029, 3.253503, 
                              3.273079, 3.451232, 3.54805, 3.6008]
iso.add_to_hdf5(fname)





