import sedonalib as sed
import numpy as np

vmin = 1e9
vmax = 2e9
ve   = 2.5e8
texp = 20.0*3600.0*24
rho0 = 1e-14
Lum  = 1e43
name = "atmosphere.mod"

mymod = sed.new_model('1D_sphere',dims=50,rmin=vmin*texp,rmax=vmax*texp)

r = mymod.rout
mymod.set_velocity(r/texp)
mymod.set_density(rho0*np.exp(-1.0*(r/texp-vmin)/ve))
mymod.add_elements(["14.28","20.40"])
mymod.set_constant_composition([0.99,0.01])
mymod.set_time(texp)

r = mymod.r
Tph = ((Lum/2/(4.0*3.1415*(vmin*texp)**2.0*5.6704e-5))**0.25)
Tph = 7700.0
T = Tph*(1 - (1 - (vmin*texp/r)**2.0)**0.5)**0.25
mymod.set_temperature(T)

mymod.plot()


mymod.write(name)
exit(0)

par = sed.new_param_file()
par.set_template("supernova/1D_typeIa_supernova_spectrum")

#par["data_atomic_file"] = "sedona_home..\"/data/cmfgen_levelcap100.hdf5\""
par["core_radius"]     = vmin*texp
par["core_luminosity"] = 1e43
par["particles_n_emit_radioactive"] = 0
par["particles_last_iter_pump"] = 1
par["core_n_emit"]     = 1e6
par["model_file"]      = name
par["opacity_bound_bound"] = 1
par["opacity_bound_free"] = 1
par["opacity_free_free"]  = 1

par["output_write_atomic_levels"] = 1

par["transport_steady_iterate"] = 2
par["opacity_use_nlte"] = 1
par["opacity_atoms_in_nlte"] = {20}
par["output_write_radiation"] = 1

par["opacity_line_expansion"] = 0
par["tstep_time_start"] = texp
par["transport_nu_grid"] = {0.8e14,5.0e16,0.001,1}
par["line_velocity_width"] = 3e10*0.005
par.write()
print(par)
