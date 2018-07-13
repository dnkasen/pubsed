import matplotlib
from matplotlib.pyplot import *
import constants as const

num_zones_bulk = 512
inner_radius = 1.e14
bulk_outer_radius = 5.e14

nH = 4.2e12 # in innermost radius
rho_power = 2.
tail_power = 10.
L = 1.e45 # for determining initial guess of temperature profile

n_elements = 5
nHe_over_nH = 0.1 # See http://nova.astro.umd.edu/Tlusty2002/solar-abun.html
nC_over_nH = 0.000331
nN_over_nH = 0.0000832
nO_over_nH = 0.000676
nSi_over_nH = 0.

fout = open("tde_static_1D_atmosphere.mod","w")


### Constants ###

A_H = 1.0
A_He = 4.0
A_C = 12.0
A_N = 14.0
A_O = 16.00
A_Si = 28.00


### Begin setup ###

dr = (bulk_outer_radius - inner_radius)/float(num_zones_bulk)

weight_normalization = 1.+ nHe_over_nH + nC_over_nH + nN_over_nH + nO_over_nH + nSi_over_nH
weight_sum = 1. * A_H + nHe_over_nH * A_He + nC_over_nH * A_C + nN_over_nH * A_N + nO_over_nH * A_O + nSi_over_nH * A_Si
mean_atomic_weight = (weight_sum)/(weight_normalization)
print "mean atomic weight is %f" % mean_atomic_weight

f_e =  (1. + nHe_over_nH * 2. + nC_over_nH * 4. + nN_over_nH * 5. + nO_over_nH * 6. + nSi_over_nH * 10.)  # This is approximate. The 6 is from the fact that, through most of the domain, oxygen will contribute 6 free electrons (two remain bound). But it doesn't really make a big difference if you use 8 instead of 6

print "electron number density over nH is %f" % f_e

tail_outer_radius = pow(2. * f_e * nH * pow(bulk_outer_radius/inner_radius,-1. * rho_power) * const.sigma_t * pow(bulk_outer_radius,tail_power)/(tail_power - 1.),1./(tail_power - 1))
#tail_outer_radius = bulk_outer_radius
num_zones_tail = int((tail_outer_radius - bulk_outer_radius)/dr)

num_bulk_radii = num_zones_bulk + 1
num_tail_radii = num_zones_tail

num_total_radii = num_bulk_radii + num_tail_radii
num_zones_total = num_zones_bulk + num_zones_tail

print "Total number of zones is ", num_zones_total

### Header ###

fout.write("1D_sphere standard\n")
fout.write(str(num_zones_total) + " " + str(inner_radius) + " " + str(0.) + " " + str(n_elements) + "\n")
#fout.write("1.1 2.4 6.12 7.14 8.16 14.28\n")
fout.write("1.1 2.4 6.12 7.14 8.16\n")

###

radius_edges = np.zeros(num_total_radii) 
radius_centers = np.zeros(num_zones_total) 
density_edges = np.zeros(num_total_radii) 
density_centers = np.zeros(num_zones_total)
temperature_edges = np.zeros(num_total_radii) 
temperature_centers = np.zeros(num_zones_total) 
thomson_tau_edges = np.zeros(num_total_radii)

T0 = pow(3. * f_e * nH * const.sigma_t * L/( (rho_power + 1.) * 4. * np.pi * const.c * inner_radius * const.a),0.25)
temperature_edges[0] = T0

print "T0 is %f" % T0

H_mass_fraction = A_H / (weight_sum) # factors of nH cancel
H_num_fraction = 1. / (weight_normalization) 
ndens0 = nH * (weight_normalization)

He_mass_fraction = A_He  * nHe_over_nH/ (weight_sum) # factors of nH cancel
C_mass_fraction = A_C * nC_over_nH / (weight_sum) # factors of nH cancel
N_mass_fraction = A_N * nN_over_nH / (weight_sum) # factors of nH cancel
O_mass_fraction = A_O * nO_over_nH / (weight_sum) # factors of nH cancel
Si_mass_fraction = A_Si * nSi_over_nH / (weight_sum) # factors of nH cancel

rho0 = const.m_p * nH * (weight_sum)
density_edges[0] = rho0

print "rho0 is %e, H mass fraction of total is %e, H number fraction of total is %e, total number density of nuclei is %e" % (rho0,H_mass_fraction,H_num_fraction, ndens0)

if (rho_power == 3):
    analytic_mass = 4. * np.pi* rho0 * pow(inner_radius,3) * log(bulk_outer_radius/inner_radius)
else: 
    analytic_mass = 4. * np.pi* rho0 * pow(inner_radius,3) * 1./(3. - rho_power) * (pow(bulk_outer_radius/inner_radius, 3. - rho_power ) - 1.)

radius_edges[0] = inner_radius
thomson_tau = 0.
bulk_mass = 0.
tail_mass = 0.

for i  in range(num_zones_bulk):
    v_out = 0.

    radius_edges[i+1] = inner_radius + (i+1) * dr
    radius_centers[i] = inner_radius + (i+0.5) * dr

    density_edges[i+1] = rho0 * pow(radius_edges[i+1]/inner_radius,-1. * rho_power)
    density_centers[i] = 0.5 * (density_edges[i] + density_edges[i+1])

    current_nH = density_centers[i]/((weight_sum) * const.m_p)

    temperature_edges[i+1] = T0 * pow(radius_edges[i+1]/inner_radius,-1. * (rho_power + 1.)/4.)
    temperature_centers[i] = 0.5 * (temperature_edges[i] + temperature_edges[i+1])
    
#    line = "%10.8e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n" % (radius_edges[i+1],v_out,density_centers[i],temperature_centers[i],H_mass_fraction,He_mass_fraction,C_mass_fraction,N_mass_fraction,O_mass_fraction,Si_mass_fraction)
    line = "%10.8e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n" % (radius_edges[i+1],v_out,density_centers[i],temperature_centers[i],H_mass_fraction,He_mass_fraction,C_mass_fraction,N_mass_fraction,O_mass_fraction)
    fout.write(line)

    bulk_mass += 4. * np.pi /3. * ( pow(radius_edges[i+1],3.) - pow(radius_edges[i],3.)) * density_centers[i] 

    thomson_tau += (radius_edges[i+1] - radius_edges[i]) * const.sigma_t * current_nH * (f_e)
    thomson_tau_edges[i+1] = thomson_tau

for i in range(num_zones_tail):
    v_out = 0.

    radius_edges[num_zones_bulk + i + 1] = bulk_outer_radius + (i+1) * dr
    radius_centers[num_zones_bulk + i] = bulk_outer_radius + (i+0.5) * dr

    density_edges[num_zones_bulk + i + 1] = density_edges[num_zones_bulk] * pow(radius_edges[num_zones_bulk + i + 1]/bulk_outer_radius,-1. * tail_power)
    density_centers[num_zones_bulk + i] = 0.5 * (density_edges[num_zones_bulk + i] + density_edges[num_zones_bulk + i+1])

    current_nH = density_centers[num_zones_bulk + i]/((1. + 4. * nHe_over_nH + 16. * nO_over_nH ) * const.m_p)

    temperature_edges[num_zones_bulk + i+1] = T0 * pow(radius_edges[num_zones_bulk + i+1]/inner_radius,-1. * (rho_power + 1.)/4.) # not changing this to reflect tail density
    temperature_centers[num_zones_bulk + i] = 0.5 * (temperature_edges[num_zones_bulk + i] + temperature_edges[num_zones_bulk + i+1])
    
#    line = "%10.8e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n" % (radius_edges[num_zones_bulk + i+1],v_out,density_centers[num_zones_bulk + i],temperature_centers[num_zones_bulk + i], H_mass_fraction,He_mass_fraction,C_mass_fraction,N_mass_fraction,O_mass_fraction,Si_mass_fraction)
    line = "%10.8e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n" % (radius_edges[num_zones_bulk + i+1],v_out,density_centers[num_zones_bulk + i],temperature_centers[num_zones_bulk + i], H_mass_fraction,He_mass_fraction,C_mass_fraction,N_mass_fraction,O_mass_fraction)
    fout.write(line)

    tail_mass += 4. * np.pi /3. * ( pow(radius_edges[num_zones_bulk + i+1],3.) - pow(radius_edges[num_zones_bulk + i],3.)) * density_centers[num_zones_bulk + i] 

    thomson_tau += (radius_edges[num_zones_bulk + i+1] - radius_edges[num_zones_bulk + i]) * const.sigma_t * current_nH * (f_e)
    thomson_tau_edges[num_zones_bulk + i+1] = thomson_tau
    


fout.close()

print "Bulk mass should be " + str(analytic_mass)
print "Bulk mass is " + str(bulk_mass) + ", which is " + str(bulk_mass/const.m_sun) + " solar"
print "Tail mass is " + str(tail_mass) + ", which is " + str(tail_mass/const.m_sun) + " solar"
print "Total mass is " + str(bulk_mass + tail_mass) + ", which is " + str((bulk_mass + tail_mass)/const.m_sun) + " solar"
print "Outer radius is " + str(radius_edges[-1])
print "Thomson optical depth is approximately " + str(thomson_tau)


loglog(radius_centers, density_centers,'-o',markersize=3)

xlabel('height (cm)')
ylabel('density (g/cm^3)')

show()
close()

thomson_tau_edges = thomson_tau * np.ones(len(thomson_tau_edges)) - thomson_tau_edges

loglog(radius_edges, thomson_tau_edges,'-o',markersize=3)

xlabel('height (cm)')
ylabel('scattering optical depth')

show()
close()


