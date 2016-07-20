import h5py
import pylab as py
import constants as pc

#############################
def calculate_milne2(temp,E,cs,gl,gc):

        nu_t = E[0]*pc.ev_to_ergs/pc.h
        nu_0 = nu_t
        cs_0 = cs[0]
        BB_0 = 2*pc.h*nu_t**3/pc.c**2
        rsum  = 0.0

        lamb = (pc.h**2/(2*pc.pi*pc.m_e*pc.k*temp))**0.5
        nc_over_ni = 2.0/lamb**3*gc/gl #*py.exp(-1.0*E[0]/pc.k_ev/temp)

       
	for i in range(1,len(E)):
	        nu_1 = E[i]*pc.ev_to_ergs/pc.h
 	        nu_m = 0.5*(nu_1 + nu_0)
    	        dnu  = nu_1 - nu_0
    	        zeta = pc.h*(nu_1 - nu_t)/pc.k/temp
    	        BB_1 = 2*pc.h*nu_1**3/pc.c**2*py.exp(-zeta)
    	        BB_m = 0.5*(BB_0 + BB_1)
    	        cs_m = 0.5*(cs[i] + cs_0)

    	        rsum += cs_m/pc.h/nu_m*dnu*BB_m
    	        nu_0 = nu_1
    	        BB_0 = BB_1
    	        cs_0 = cs[i]
        rsum *= 4.0*pc.pi/nc_over_ni
	print temp,E[0],lamb,rsum,nc_over_ni
 	return rsum

#############################
def calculate_milne(temp,E,cs,gl,gc):

  # Maxwell-Bolztmann constants
  v_MB = (2*pc.k*temp/pc.m_e)**0.5
  MB_A = 4.0/(pc.pi)**0.5*v_MB**(-3.0)
  MB_B = pc.m_e/pc.k/2.0/temp
  milne_fac = (pc.h/pc.c/pc.m_e)**2.0

  # starting values
  sum   = 0.0
  nu_t  = E[0]*pc.ev_to_ergs/pc.h
  nu    = nu_t;
  vel   = 0.0
  fMB   = 0.0
  sigma = 0.0
  coef  = 0.0
  old_vel  = vel;
  old_coef = coef

#   // integrate over velocity/frequency
  for i in range(1,len(E)):
    nu  = E[i]*pc.ev_to_ergs/pc.h
    vel = (2*pc.h*(nu - nu_t)/pc.m_e)**0.5
    fMB   = MB_A*vel*vel*py.exp(-MB_B*vel*vel)
    sigma = milne_fac*cs[i]*nu*nu/vel/vel
    coef  = vel*sigma*fMB
    if (nu < nu_t): coef = 0.0

#     // integrate
    sum += 0.5*(coef + old_coef)*(vel - old_vel)

#     // store old values
    old_vel  = vel
    old_coef = coef
  
  return (1.0*gl)/(1.0*gc)*sum

#############################
def calculate_hydrogenic_cs(Ex):
	cs = 6.0e-18*(Ex/Ex[0])**(-3.0)
	return cs			

#############################
## 3 level plus continum hydrogen
#############################

fname = 'H_3lev_plusC_atomdata.hdf5'
level_g = [2,    8,   18,   1]
level_E = [0, 10.2, 12.09,  0]
level_i = [0,    0,    0,   1]
lu = [  1,  2,  2]
ll = [  0,  0,  1]
Aij = [4.696e8, 5.572e7, 4.408e7] 
ion_chi    = [13.6, 999999]
ion_ground = [0, 3]

#############################
## 2 level no real continum hydrogen
#############################

#fname = 'H_2lev_noC_atomdata.hdf5'
#level_g = [2,    8]
#level_E = [0, 10.2]
#level_i = [0,    0]
#lu = [  1]
#ll = [  0]
#Aij = [1.0e9]
#ion_chi    = [13.6]
#ion_ground = [0]


#############################
## 2 ion stages
#############################

#fname = 'H_1lev_2ion_atomdata.hdf5'
#level_g = [2,    1]
#level_E = [0,    0]
#level_i = [0,    1]
#lu = [  ]
#ll = [  ]
#Aij = []
#ion_chi    = [13.6,99999]
#ion_ground = [0,1]

#############################
#############################
#############################
## write data file

# grid for recombination coefficients
t_start = 1.0e2
t_stop  = 1.0e8
delta_T = 0.1


f = h5py.File(fname,'w')

base = "1/"
grp = f.create_group(base)

f[base].attrs["n_ions"]    = len(ion_chi)
f[base].attrs["n_levels" ] = int(len(level_g))
f[base].attrs["n_lines" ]  = len(Aij)

f.create_dataset(base + "ion_chi",data = ion_chi)
f.create_dataset(base + "ion_ground",data = ion_ground)
f.create_dataset(base + "level_g",data = level_g)
f.create_dataset(base + "level_E",data = level_E)
f.create_dataset(base + "level_i",data = level_i)
f.create_dataset(base + "line_l",data = ll)
f.create_dataset(base + "line_u",data = lu)
f.create_dataset(base + "line_A",data = Aij)

photo = base +'/bound_free/'
f.create_group(photo)

# temperature array
t_arr = []
t = t_start
while (t <= t_stop):
	t_arr.append(t)
	t = t + t*delta_T
t = t + t*delta_T
t_arr.append(t)

milne = py.zeros(len(t_arr))


npts = 10000
fmax = 5
for i in range(len(level_E)):

	# ionization potential
	ion  = level_i[i]
	chi  = ion_chi[ion]

	# energy range
	Eion = chi - level_E[i]
	Emax = fmax*Eion
	dE = (Emax - Eion)/(1.0*npts)
	Ex = py.arange(Eion,Emax,dE)

	# cross-sections
	cs = calculate_hydrogenic_cs(Ex)

	# statistical weights
	gl = level_g[i]
	ion_to = ion+1
	if (ion_to > len(ion_chi)): gc = 1
	else: gc = level_g[ion_ground[ion_to]]
	
        calculate_milne2(50000.0,Ex,cs,gl,gc)
        exit(0)

	# milne
	for j in range(len(t_arr)):
		milne[j] = calculate_milne2(t_arr[j],Ex,cs,gl,gc)

	tref,aref = py.loadtxt('alpha_1.dat',unpack=1)
	py.plot(tref,aref,color='red')

	py.plot(t_arr,milne)
	py.yscale('log')
	py.xscale('log')
	py.ion()
	py.show()
	j = raw_input()

	# write out
	f[base].attrs[photo + str(i) + "_photoi_npts"]    = len(Ex)
	f.create_dataset(photo + str(i) + '_photoi_E',data = Ex)
	f.create_dataset(photo + str(i) + '_photoi_s',data = cs)
#	f[base].attrs[photo + str(i) + "_recomb_npts"]    = len(t_arr)	
#	f.create_dataset(photo + str(i) + '_recomb_a',data = milne)
#	f.create_dataset(photo + str(i) + '_recomb_T',data = t_arr)



f.close()

