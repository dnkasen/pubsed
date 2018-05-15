'''

This script is used to write the hdf5 atomic data files used by sedona6.
It reads the raw cmfgen data from "cmfgen_data.hdf5" and writes another
hdf5 file that contains only the species contained in the config.yml file.

It has various options for processing the raw data. You can control these
by editing the "merge_denegerate", "skip_autoionizing" and "level_cap"
variables below. 

To control which atoms / species go in to your model, edit the config.yml
file. 

Usage: python create_model.py

'''

import yaml
import h5py
from atom import Atom


# merge levels with exactly the same energy
merge_degenerate = True

# skip levels with energies greater than the ionization energy for the species
skip_autoionizing = True

# use at most this many levels per species.
level_cap = 100

# which config file to use
config_file = "config.yml"


with open(config_file, 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

with h5py.File('cmfgen_data.hdf5', 'r') as h5f:
    atoms = []
    for atom in cfg:
        print(atom)
        atoms.append(Atom(h5f, cfg[atom], merge_degenerate, skip_autoionizing, level_cap))

with h5py.File('atom_data.hdf5', 'w') as h5f:
    for atom in atoms:
        atom.write_hdf5(h5f)
