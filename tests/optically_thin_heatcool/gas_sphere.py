#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import constants as const
import numpy as np
from math import pow, sqrt
import blackbody_emitter as bb
import element


class GasSphere:

#params (default values)
    def __init__(self):

        self.min_radius = 3.e13
        self.max_radius = 1.e15
        self.density_power = 0.
        self.peak_number_density = 1.e1

        self.num_radius_centers = 128
        self.num_radius_edges = self.num_radius_centers + 1
    
        self.radius_edges = []
        self.radius_edges_nondim = []

        self.radius_centers = []
        self.radius_centers_nondim = []

        self.elements = []
        self.element_number_fractions = []

#################
    
    def set_num_radii(self, nr):
        self.num_radius_centers = nr
        self.num_radius_edges = nr + 1

################

# Add an element with hydrogenic approximations for the effective quantum "n" for each level, which are used for  the hydrogenic photoionization cross sections for each level. Future work will need to go beyond these approximations.

# You'll also need to keep track of ionization states for elements other than hydrogen
    def add_element(self,Z,number_fraction,level_energy_list,level_g_list):
        this_element = element.Element(Z)

        for level_index in range(len(level_energy_list)):

            n_eff = pow(this_element.total_ionization_energy_ev/(this_element.total_ionization_energy_ev - level_energy_list[level_index]),0.5) # this needs to be modified depending on the ionizaiton state for elements other than hydrogen
            this_element.add_level(level_index,n_eff,level_energy_list[level_index],level_g_list[level_index])
        
        self.elements.append(this_element)




    def n_r(self,r): # Total number density of nuclei (I'm imagining this includes sum over all elements)
        return self.peak_number_density * pow(r/self.min_radius,-1. * self.density_power)


    def create_radial_grids(self): 
        self.num_radius_edges = self.num_radius_centers+1
        self.create_radial_edge_grid()
        self.create_radial_center_grid()

    def create_radial_edge_grid(self):
        self.radius_edges = np.linspace(self.min_radius,self.max_radius,self.num_radius_edges)

    def create_radial_center_grid(self):
        self.radius_centers = np.fromiter( (0.5 * (self.radius_edges[index] + self.radius_edges[index+1] ) for index in range(self.num_radius_centers)), np.float)



#    def nH_nondim(self, r):
#        return pow(r/self.min_radius,-1. * self.density_power)

#    def nHe_nondim(self, r):
#        return self.n_H(r) * self.helium_to_hydrogen_number_fraction 

# oxygen ...

#################
#    def IonizationDepletion(self, r):
#        return atoms.alpha_Bs[0] * pow(self.nH_nondim(r) * self.peak_number_density * r,2)





