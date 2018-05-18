import os
import yaml
import h5py
import numpy as np
from species import Species

inv_cm_to_ev = 0.00012398424

class Atom:

    def __init__(self, h5f, data_dict, merge_degenerate=True,
                 skip_autoionizing=True, level_cap=None):

        self.Z = data_dict['Z']

        self.species = []
        for spec in data_dict['species']:
            group = os.path.join(str(self.Z), str(spec))
            self.species.append(Species(h5f, group, merge_degenerate,
                                        skip_autoionizing, level_cap))

        # ionization energies - should convert these to eV
        ion_chi = [s.ion_chi for s in self.species]
        ion_chi += [1.0e8 / inv_cm_to_ev]
        self.ion_chi = np.array(ion_chi)
        
        # the ground state of each ionic species
        ion_ground = [0] + [s.nlevels for s in self.species]
        ion_ground = np.array(ion_ground)
        self.ion_ground = np.cumsum(ion_ground)

        # the ionic species of each level
        level_i = np.concatenate([i*np.ones_like(s.E) for i, s in enumerate(self.species)])
        level_i = np.concatenate([level_i, [len(self.species)]])
        self.level_i = np.array(level_i)
        
        # the energies of each level
        level_E = np.concatenate([s.E for s in self.species])
        level_E = np.concatenate([level_E, [0.0]])
        self.level_E = level_E

        # the degeneracies of each level
        level_g = np.concatenate([s.g for s in self.species])
        level_g = np.concatenate([level_g, [0.0]])
        self.level_g = level_g

        line_l = np.concatenate([off + s.line_l for off, s in zip(self.ion_ground, self.species)])
        line_u = np.concatenate([off + s.line_u for off, s in zip(self.ion_ground, self.species)])
        line_A = np.concatenate([s.line_A for s in self.species])

        self.line_l = line_l
        self.line_u = line_u
        self.line_A = line_A

        self.level_E *= inv_cm_to_ev
        self.ion_chi *= inv_cm_to_ev

        assert(np.all(self.line_l < self.line_u))

    def write_hdf5(self, h5f):

        grp = h5f.create_group(str(self.Z))

        grp.create_dataset("ion_chi", data=self.ion_chi)
        grp.create_dataset("ion_ground", data=self.ion_ground)

        grp.create_dataset("level_i", data=self.level_i)
        grp.create_dataset("level_E", data=self.level_E)
        grp.create_dataset("level_g", data=self.level_g)

        grp.create_dataset("line_l", data=self.line_l)
        grp.create_dataset("line_u", data=self.line_u)
        grp.create_dataset("line_A", data=self.line_A)

        grp.attrs['n_ions']   = len(self.ion_chi)
        grp.attrs['n_levels'] = len(self.level_E)
        grp.attrs['n_lines']  = len(self.line_A)
        

if __name__ == "__main__":

    config_file = "config.yml"
    with open(config_file, 'r') as ymlfile:
        cfg = yaml.load(ymlfile)

    with h5py.File('cmfgen_data.hdf5', 'r') as h5f:
        atoms = []
        for atom in cfg:
            print(atom)
            atoms.append(Atom(h5f, cfg[atom]))

    with h5py.File('atom_data.hdf5', 'w') as h5f:
        for atom in atoms:
            atom.write_hdf5(h5f)
