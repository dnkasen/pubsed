import h5py
import numpy as np
from collections import defaultdict

from util import to_array, iteritems


class Species:

    def __init__(self, h5f, species_string, 
                 merge_degenerate=True,
                 skip_autoionizing=True,
                 level_cap=None):
        '''

        This is class for representing the atomic line data for a
        single Species (e.g. CII).

        '''

        self.merge_degenerate_levels = merge_degenerate
        self.skip_autoionizing_levels = skip_autoionizing
        self.level_cap = level_cap
        
        self._read_data(h5f, species_string)

        if self.skip_autoionizing_levels or self.level_cap:
            self._cap_levels()

        if self.merge_degenerate_levels:
            self._merge_degenerate()

    def _read_data(self, h5f, species_string):

        self.nlevels = h5f[species_string].attrs['n_levels']
        self.nlines  = h5f[species_string].attrs['n_lines']
        self.ion_chi = h5f[species_string].attrs['ion_chi']
        
        self.E = h5f[species_string]['level_E'][:]
        self.g = h5f[species_string]['level_g'][:]
        self.lev = h5f[species_string]['level_i'][:]
        
        self.max_level = self.nlevels - 1
        
        self.line_l = h5f[species_string]['line_l'][:]
        self.line_u = h5f[species_string]['line_u'][:]
        self.line_A = h5f[species_string]['line_A'][:]

    def _cap_levels(self):
        '''

        This routine removes the highest energy levels from the atom, based on
        the values of self.max_level and self.skip_autoionizing.


        '''

        if self.skip_autoionizing_levels:
            first_bad_index = np.argmax(self.E > self.ion_chi)
        else:
            first_bad_index = len(self.E)
        
        if self.level_cap:
            if first_bad_index == 0:
                first_bad_index = self.level_cap
            else:
                first_bad_index = min(first_bad_index, self.level_cap)

        if (first_bad_index == 0):
            return

        self.E = self.E[:first_bad_index]
        self.g = self.g[:first_bad_index]
        self.lev = self.lev[:first_bad_index]
        
        locs = self.line_u < first_bad_index
        self.line_l = self.line_l[locs]
        self.line_u = self.line_u[locs]
        self.line_A = self.line_A[locs]

        self.nlevels = len(self.E)
        self.max_level = self.nlevels - 1

    def _merge_degenerate(self):
        '''

        This removes degenerate energy levels from the atomic data by merging them
        into one superlevel. We treat this is a graph problem, treating the energy
        levels as vertices and the transitions as edges. Merging two energy levels
        then translates into performing vertex contraction on the graph.

        '''

        # convert from arrays to adjacency list
        graph = defaultdict(set)
        weights = defaultdict(float)
        for l, u, A in zip(self.line_l, self.line_u, self.line_A):
            graph[l].add(u)
            graph[u].add(l)
            weights[frozenset([l, u])] = A

        def merge_levels(l, u):
            '''
            
            This function merges level u into level l, averaging the energy and summing 
            the degeneracy factors. The Einstein A coefficients will be combined using
            a weighted average based on the degeneracy factors.
            
            '''

            g_tot =      self.g[l-1] + self.g[u-1]
            E_avg = 0.5*(self.E[l-1] + self.E[u-1])
            
            # for all edges going out of level l, replace A with a 
            # degeneracy factor weighted average over levels u and l.
            for v in graph[l]:
                if (v == u): continue
                key_l = frozenset([l, v])
                key_u = frozenset([u, v])
                weights[key_l] = (self.g[l-1]*weights[key_l] + 
                                  self.g[u-1]*weights[key_u])/g_tot
                
            # replace level l's g and E with the combined values
            self.g[l-1] = g_tot
            self.E[l-1] = E_avg
        
            # give level l the combination of all edges going into levels l and u
            combined_vertices = graph[l] | graph[u]
            combined_vertices.discard(u)
            combined_vertices.discard(l)
            graph[l] = combined_vertices
            
            # delete level u from the graph
            graph.pop(u)

            # remove u and add l to all neighbors
            for v in combined_vertices:
                graph[v].discard(u)
                graph[v].add(l)

            assert(u not in graph)
            for k, v in iteritems(graph):
                assert(u not in v)

        # now perform the merge for every set of degenerate levels
        def find_degenerate_levels(E, lev):
            _, indices = np.unique(E, return_index=True)
            indices = np.sort(indices)  # surprised this isn't always true
            unique_levs = np.split(lev, indices[1:])
            degenerate_levels = np.array([levs for levs in unique_levs if levs.size > 1])
            return degenerate_levels, indices

        degenerate_levels, indices = find_degenerate_levels(self.E, self.lev)

        if degenerate_levels.size == 0:  # no levels will be merged
            return

        for dlev in degenerate_levels:
            assert(len(dlev) > 1)
            l = dlev[0]
            for u in dlev[1:]:
                merge_levels(l, u)

        # convert back to array representation of line transitions
        new_l, new_u, new_A = [], [], []
        for k, vals in iteritems(graph):
            for v in vals:
                assert(k != v)
                if (k < v):
                    new_l.append(k)
                    new_u.append(v)
                    new_A.append(weights[frozenset([k, v])])

        new_l, new_u, new_A = to_array(new_l, new_u, new_A)
        new_lev = self.lev[indices]
        
        # now re-number the levels so that they are continuous from 0 to nlevels
        offset = indices - np.arange(indices.size)
        mapping = np.digitize(self.lev, indices, right=True) 
        
        new_l   -= offset[mapping[new_l]]
        new_u   -= offset[mapping[new_u]]
        new_lev -= offset[mapping[new_lev]]
        assert(np.all(new_lev == np.arange(new_lev.size)))

        # go ahead and sort too
        new_l, new_u, new_A = (np.array(t) for t in zip(*sorted(zip(new_l, new_u, new_A))))

        self.line_l = new_l
        self.line_u = new_u
        self.line_A = new_A

        self.lev = new_lev
        self.E = self.E[indices]
        self.g = self.g[indices]

        self.nlevels = len(self.lev)
        self.max_level = self.nlevels - 1
        self.nlines = len(self.line_l)

        assert(np.all(self.line_l < self.line_u))


if __name__ == "__main__":

    data_file = "cmfgen_data.hdf5"

    with h5py.File(data_file, "r") as h5f:
        species = Species(h5f, '7/2')  # singly ionized nitrogen
    
        print(species.E)
        print(species.g)
        print(species.lev)

        print(species.line_l)
        print(species.line_u)
        print(species.line_A)
        
