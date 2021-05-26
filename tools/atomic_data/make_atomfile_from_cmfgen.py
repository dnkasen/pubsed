'''

This script is used to convert the raw cmfgen data, which comes in
various ascii formats, into a standardized hdf5 file.

To configure this, you can edit the 'base_dir' and 'data_sources'
variables below.

    Usage: python make_atomfile_from_cmfgen.py

'''

import os
import re
import numpy as np
import h5py

from util import \
    peek_line, \
    skip_until_found, \
    to_array, \
    iteritems, \
    open_text_file


# this should point to the directory containing the raw cmfgen data
# on your machine, (i.e. the directory containing the ascii files)
base_dir = '/Users/kasen/data/cmfgen_atomic_raw/'


# this dictionary specifies which data files to use for each species.
# Usually, it is recommended to use the most recent files, but for
# some species, like Hydrogen, the older data is recommended. See
# See the README files in the species sub-directories in CMFGEN
# for more information.
data_sources = {
     '13':
           {7 : 'AL/VII/23oct02/fin_osc',
            11: 'AL/XI/23oct02/fin_osc',
            8:  'AL/VIII/23oct02/fin_osc',
            6:  'AL/VI/23oct02/fin_osc',
            2:  'AL/II/5aug97/al2_osc.dat',
            9:  'AL/IX/23oct02/fin_osc',
            1:  'AL/I/29jul10/fin_osc',
            4:  'AL/IV/23oct02/fin_osc',
            10: 'AL/X/23oct02/fin_osc',
            3:  'AL/III/30oct12/aliii_osc_split.dat',
            5:  'AL/V/23oct02/fin_osc'},
     '18':
           {7: 'ARG/VII/15feb01/arviiosc_rev.dat',
            8: 'ARG/VIII/15feb01/arviiiosc_rev.dat',
            6: 'ARG/VI/15feb01/arviosc_rev.dat',
            2: 'ARG/II/9sep11/fin_osc',
            1: 'ARG/I/9sep11/fin_osc',
            4: 'ARG/IV/1dec99/fin_osc.dat',
            3: 'ARG/III/19nov07/fin_osc',
            5: 'ARG/V/1dec99/arvosc_rev.dat'},
     '20':
           {7:  'CA/VII/10apr99/osc_op_sp.dat',
            8:  'CA/VIII/11jun01/caviii_osc.dat',
            6:  'CA/VI/10apr99/osc_op_sp.dat',
            2:  'CA/II/30oct12/ca2_osc_split.dat',
            9:  'CA/IX/11jun01/caix_osc.dat',
            1:  'CA/I/5aug97/cai_osc_split.dat',
            4:  'CA/IV/10apr99/osc_op_sp.dat',
            10: 'CA/X/30mar02/cax_osc.dat',
            3:  'CA/III/10apr99/osc_op_sp.dat',
            5:  'CA/V/10apr99/osc_op_sp.dat'},
     '6':
           {6: 'CARB/VI/9may02/fin_osc',
            2: 'CARB/II/30oct12/c2osc_rev.dat',
            1: 'CARB/I/12dec04/ci_split_osc',
            4: 'CARB/IV/30oct12/civosc_a12.dat',
            3: 'CARB/III/23dec04/ciiiosc_st_split_big.dat',
            5: 'CARB/V/20oct02/fin_osc_split'},
     '17':
           {7: 'CHL/VII/15feb01/clviiosc_rev.dat',
            6: 'CHL/VI/15feb01/clviosc_rev.dat',
            4: 'CHL/IV/15feb01/clivosc_fin.dat',
            5: 'CHL/V/15feb01/clvosc_fin.dat'},
     '24':
           {6: 'CHRO/VI/18oct00/crvi_osc.dat',
            2: 'CHRO/II/15aug12/crii_osc.dat',
            1: 'CHRO/I/10aug12/cri_osc.dat',
            4: 'CHRO/IV/18oct00/criv_osc.dat',
            3: 'CHRO/III/18oct00/criii_osc.dat',
            5: 'CHRO/V/18oct00/crv_osc.dat'},
     '27':
           {7: 'COB/VII/18oct00/covii_osc.dat',
            8: 'COB/VIII/18oct00/coviii_osc.dat',
            6: 'COB/VI/18oct00/covi_osc.dat',
            2: 'COB/II/15nov11/fin_osc_bound',
            9: 'COB/IX/18oct00/coix_osc.dat',
            4: 'COB/IV/4jan12/coiv_osc.dat',
            3: 'COB/III/30oct12/coiii_osc.dat',
            5: 'COB/V/18oct00/cov_osc.dat'},
     '26':
           {7:  'FE/VII/18oct00/fevii_osc.dat',
            11: 'FE/XI/16oct02/osc_split_op',
            8:  'FE/VIII/8may97/feviii_osc_kb_rk.dat',
            6:  'FE/VI/18oct00/fevi_osc.dat',
            14: 'FE/XIV/14jan08/osc_split_op',
            2:  'FE/II/10sep16/fe2_osc',
            9:  'FE/IX/25aug03/feix_osc',
            13: 'FE/XIII/14jan08/osc_split_op',
            12: 'FE/XII/16oct02/osc_split_op',
            1:  'FE/I/07sep16/fei_osc',
            15: 'FE/XV/14jan08/osc_split_op',
            4:  'FE/IV/18oct00/feiv_osc_rev2.dat',
            16: 'FE/XVI/14jan08/osc_split_op',
            10: 'FE/X/5may02/osc_split_op',
            3:  'FE/III/30oct12/FeIII_OSC',
            5:  'FE/V/18oct00/fev_osc.dat'},
     '2':
           {1: 'HE/I/15jul15/hei_osc',
            2: 'HE/II/5dec96/he2_osc.dat'},
     '1':
           {1: 'HYD/I/5dec96/hi_osc.dat'},
     '25':
           {7: 'MAN/VII/18oct00/mnvii_osc.dat',
            6: 'MAN/VI/18oct00/mnvi_osc.dat',
            2: 'MAN/II/18oct00/mnii_osc.dat',
            4: 'MAN/IV/18oct00/mniv_osc.dat',
            3: 'MAN/III/18oct00/mniii_osc.dat',
            5: 'MAN/V/18oct00/mnv_osc.dat'},
     '12':
           {7:  'MG/VII/20jun01/mgviiosc_rev.dat',
            8:  'MG/VIII/23oct02/fin_osc',
            6:  'MG/VI/20jun01/mgviosc_rev.dat',
            2:  'MG/II/30oct12/mg2_osc_split.dat',
            9:  'MG/IX/23oct02/fin_osc',
            1:  'MG/I/5aug97/mgi_osc_split.dat',
            4:  'MG/IV/20jun01/mgivosc_rev.dat',
            10: 'MG/X/23oct02/fin_osc',
            3:  'MG/III/20jun01/mgiiiosc_rev.dat',
            5:  'MG/V/20jun01/mgvosc_rev.dat'},
     '11':
           {7: 'NA/VII/20jun01/naviiosc_rev.dat',
            8: 'NA/VIII/23oct02/fin_osc',
            6: 'NA/VI/20jun01/naviosc_rev.dat',
            2: 'NA/II/15feb01/na2_osc.dat',
            9: 'NA/IX/23oct02/fin_osc',
            1: 'NA/I/5aug97/nai_osc_split.dat',
            4: 'NA/IV/15feb01/naivosc_rev.dat',
            3: 'NA/III/15feb01/naiiiosc_rev.dat',
            5: 'NA/V/20jun01/navosc_rev.dat'},
     '10':
           {7: 'NEON/VII/19nov04/neviiosc_rev.dat',
            8: 'NEON/VIII/20jun01/neviiiosc_rev.dat',
            6: 'NEON/VI/20jun01/neviosc_rev.dat',
            2: 'NEON/II/19nov07/fin_osc',
            1: 'NEON/I/9sep11/fin_osc',
            4: 'NEON/IV/1dec99/fin_osc.dat',
            3: 'NEON/III/19nov07/fin_osc',
            5: 'NEON/V/20jun01/nevosc_rev.dat'},
     '28':
           {7: 'NICK/VII/18oct00/nkvii_osc.dat',
            8: 'NICK/VIII/11jun01/nkviii_osc.dat',
            6: 'NICK/VI/18oct00/nkvi_osc.dat',
            2: 'NICK/II/30oct12/nkii_osc.dat',
            9: 'NICK/IX/11jun01/nkix_osc.dat',
            4: 'NICK/IV/18oct00/nkiv_osc.dat',
            3: 'NICK/III/27aug12/nkiii_osc.dat',
            5: 'NICK/V/18oct00/nkv_osc.dat'},
     '7':
          {7: 'NIT/VII/4feb05/fin_osc',
           6: 'NIT/VI/4feb05/fin_osc_pack',
           2: 'NIT/II/23jan06/fin_osc',
           1: 'NIT/I/12sep12/ni_osc',
           4: 'NIT/IV/4nov10/nivosc',
           3: 'NIT/III/24mar07/niiiosc_rev.dat',
           5: 'NIT/V/30oct12/nvosc_a12_split.dat'},
     '8':
          {7: 'OXY/VII/20sep02/fin_osc_pack',
           8: 'OXY/VIII/20sep02/fin_osc_split',
           6: 'OXY/VI/30oct12/osixosc_a12_split.dat',
           2: 'OXY/II/23mar05/o2osc_fin.dat',
           1: 'OXY/I/20sep11/oi_osc_mchf',
           4: 'OXY/IV/19nov07/fin_osc',
           3: 'OXY/III/15mar08/oiiiosc',
           5: 'OXY/V/5dec96/ovosc_ns_split.dat'},
     '15':
          {2: 'PHOS/II/7oct15/p2_osc',
           4: 'PHOS/IV/15feb01/pivosc_rev.dat',
           3: 'PHOS/III/7oct15/piii_osc',
           5: 'PHOS/V/15feb01/pvosc_rev.dat'},
     '19':
          {6: 'POT/VI/15feb01/kviosc_rev.dat',
           2: 'POT/II/4mar12/fin_osc',
           1: 'POT/I/4mar12/fin_osc',
           4: 'POT/IV/15feb01/fin_osc.dat',
           3: 'POT/III/15feb01/kiiiosc_rev.dat',
           5: 'POT/V/15feb01/fin_osc.dat'},
     '21':
          {2: 'SCAN/II/01jul13/fin_osc',
           3: 'SCAN/III/3dec12/fin_osc'},
     '14':
          {7: 'SIL/VII/2dec15/sksev_osc',
           6: 'SIL/VI/30mar02/skvi_osc.dat',
           2: 'SIL/II/16sep15/si2_osc_nist',
           1: 'SIL/I/23nov11/SiI_OSC',
           4: 'SIL/IV/30oct12/osc_op_split.dat',
           3: 'SIL/III/23aug97/osc_op.dat',
           5: 'SIL/V/30mar02/skv_osc.dat'},
     '16':
          {6: 'SUL/VI/26sep14/sviosc_fin.dat',
           2: 'SUL/II/26sep14/s2_osc',
           1: 'SUL/I/26sep14/SI_OSC',
           4: 'SUL/IV/26sep14/sivosc_fin',
           3: 'SUL/III/26sep14/siiiosc_fin',
           5: 'SUL/V/26sep14/svosc_fin.dat'},
     '22':
          {2: 'TIT/II/5dec15/tkii_osc.dat',
           3: 'TIT/III/7dec15/tkiii_osc.dat',
           4: 'TIT/IV/16oct15/tkiv_osc.dat'},
}

class CollisionalData:

    def __init__(self,id):

        self.id = id
        self.type = ""
        self.O = []
        self.T = []
        self.lu = 0
        self.ll = 0

class CrossSection:

    def __init__(self,id):

        self.id = id
        self.type = "OP"
        self.E = []
        self.s = []
        self.configuration = ""


class CMFGENDataReader:

    # used to pick the line with defining the levels involved in
    # each transition out of the ascii file
    transition_matcher = re.compile("(.*?)([0-9]+-\s*[0-9]+)(.*)")

    def __init__(self, base,fn,species,ion):
        print(base + '/' + fn, species, ion)
        self.data_source = np.string_("CMFGEN: " + fn)
        self._read_data(base + '/' + fn)
        self._convert_index()
        self._read_photoion_data(species,ion)
        self._read_collisional_data("cmfgen_col_data.h5",species,ion)

    ## Read in level, line data of ion
    def _read_data(self, fn):

        inv_cm_to_ev = 0.00012398424

        assert(os.path.isfile(fn))

        with open_text_file(fn, "r") as f:
            header = self._read_header(f)

            self.nlevels = int(header['Number of energy levels'])
            self.nlines  = int(header['Number of transitions'])
            self.ion_chi = float(header['Ionization energy'])
            # convert from cm^(-1) to eV
            self.ion_chi *= inv_cm_to_ev

            assert(self.nlevels > 0)
            assert(self.nlines  >= 0)
            assert(np.all(self.ion_chi > 0))

            self.E, self.g, self.lev = [], [], []

            self.measured = []
            self.orbitals = []
            self.parity   = []

            id_index = -1
            line = f.readline()
            while line.strip() != '':

                entries = re.split('\s+', line.strip())

                orbital = entries[0]
                parity = self._get_parity(orbital)
                self.orbitals.append(np.string_(orbital))
                self.parity.append(parity)

                if id_index < 0:
                    id_index = self._get_id_index(entries)

                degeneracy = float(entries[1])
                energy = float(entries[2])

                # negative id means that the level doesn't have a measured energy
                identifier = int(entries[id_index])
                level = abs(identifier)
                self.measured.append(identifier > 0)

                self.g.append(degeneracy)
                self.E.append(energy*inv_cm_to_ev)
                self.lev.append(level)

                line = f.readline()


            self.nlevels = len(self.E)
            self.max_level = self.nlevels  # note, 1-based indexing

            self.E, self.g, self.lev = to_array(self.E, self.g, self.lev)
            self.measured, self.orbitals, self.parity = to_array(self.measured, self.orbitals, self.parity)


            # Read line transition
            if (self.nlines == 0):
                return

            skip_until_found('Transition', f)

            # Skip some blank lines
            while line.strip() == '':
                line = f.readline()

            self.line_l, self.line_u, self.line_A = [], [], []

            while line.strip() != '':
                l, u = self._get_transition_levels(line)
                A = self._get_einstein_A(line)

                if (u <= self.max_level):
                    self.line_l.append(l)
                    self.line_u.append(u)
                    self.line_A.append(A)

                line = f.readline()

            assert(len(self.line_l) == len(self.line_u) == len(self.line_A))
            self.nlines = len(self.line_l)

            self.line_l, self.line_u, self.line_A = to_array(self.line_l, self.line_u, self.line_A)




    def _read_header(self, f):
        skip_until_found('!Date', f)  # we assume there will always be a "Date"
        line = peek_line(f)
        header = {}
        while line.strip() != '':
            val, key = line.split('!')
            val, key = val.strip(), key.strip()
            header[key] = val
            line = f.readline()
        return header

    def _get_transition_levels(self, line):
        m = self.transition_matcher.match(line)
        transition_entry = m.groups()[1]
        vals = transition_entry.split('-')
        return int(vals[0]), int(vals[1])

    def _get_einstein_A(self, line):

        def is_float(s):
            try:
                _ = float(s)
                return True
            except ValueError:
                return False

        entries = re.split('\s+', line.strip())
        entries = [entry for entry in entries if is_float(entry)]
        return float(entries[1])  # always the second floating point entry

    def _get_parity(self, orbital):
        symbol = str(orbital).split("_")[-1]
        odd    = symbol.find('o') >= 0
        even   = symbol.find('e') >= 0
        if (even and odd):
            return -1
        if ((not even) and (not odd)):
            return -1
        return int(odd)

    def _get_id_index(self, entries):

        def is_integer(s):
            try:
                _ = int(s)
                return True
            except ValueError:
                return False

        results = [is_integer(s) for s in entries]

        assert(sum(results) == 1)  # should only have integer per line
        return results.index(True)

    def _convert_index(self):
        '''

        Convert from 1- to 0- based indexing.

        '''

        self.lev -= 1
        if (self.nlines > 0):
            self.line_l -= 1
            self.line_u -= 1


    ####################################
    ## Read in photoionizaiton cross-section
    ## data from file
    ####################################
    def _read_photoion_data(self,species,ion):

        self.cs_data  = []

        # default is no detailed photoionization cross section data
        self.level_cs = 0*np.arange(self.nlevels)  - 1


        # right now we only have detailed data for hydrogen
        if (int(species) == 1):
            fin = h5py.File("topbase_photo_cs/topbase_PI.hdf5","r")
            n_cs = np.array(fin["1/0/n_cs"])
            for i in range(n_cs):
                base = str(species) + "/" + str(ion) + "/cs_" + str(i) + "/"
                cs = CrossSection(i)
                cs.E = np.array(fin[base + "E_ev"])
                cs.s = np.array(fin[base + "sigma"])
                cs.type = "Topbase"
                cs.configuration = "4D(3)" + str(i)
                self.cs_data.append(cs)

                # here need to figure out which level is which cs
                self.level_cs = 0*np.arange(self.nlevels)  - 1 + n_cs
                for i in range(n_cs):
                    self.level_cs[i] = i


    ####################################
    ## Read in photoionizaiton cross-section
    ## data from file
    ####################################
    def _read_collisional_data(self,fname,species,ion):

        self.col_data  = []

        print("col data for " + str(species) + " " + str(ion))
        # open hdf5 file
        fin = h5py.File(fname,'r')
        base = str(species) + '/' + str(ion) + '/'
        if (not base + "n_data" in fin):
            return

        n_data = int(np.array(fin[base + 'n_data']))


        # faked up for now
        for i in range(n_data):
            name = base + '/' + str(i) + '/'
            col = CollisionalData(i)
            col.T = np.array(fin[name + 'T'])
            col.O = np.array(fin[name + 'C'])
            col.type = "cmfgen"
            col.ll = int(np.array(fin[name + 'lev_l']))
            col.lu = int(np.array(fin[name + 'lev_u']))
            self.col_data.append(col)

    ####################################
    ## Write ion data to hdf5 file
    ####################################
    def write_to_file(self,fname,species,ion):

        h5f = h5py.File(fname, 'a')

        if (not str(species) in h5f):
            species_group = h5f.create_group(str(species))
        ion_group = h5f.create_group(str(species) +'/' + str(ion))

        ion_group.attrs['n_levels'] = int(self.nlevels)
        ion_group.attrs['n_lines']  = self.nlines

        ion_group.create_dataset('data_source',data =self.data_source)
        ion_group.create_dataset('n_levels',data=int(self.nlevels))
        ion_group.create_dataset('n_lines',data=int(self.nlines))

        #ion_group.cattrs['ion_chi']  = self.ion_chi
        ion_group.create_dataset('ion_chi',data=self.ion_chi)

        ion_group.create_dataset("level_i", data=self.lev)
        ion_group.create_dataset("level_g", data=self.g)
        ion_group.create_dataset("level_E", data=self.E)
        ion_group.create_dataset("level_cs", data=self.level_cs)

        ion_group.create_dataset("level_measured", data=self.measured)
        ion_group.create_dataset("level_config", data=self.orbitals)
        ion_group.create_dataset("level_parity",   data=self.parity)

        if (self.nlines > 0):
            ion_group.create_dataset("line_l", data=self.line_l)
            ion_group.create_dataset("line_u", data=self.line_u)
            ion_group.create_dataset("line_A", data=self.line_A)

        # write photoionization data
        cs_group = ion_group.create_group("photoion_data")
        cs_group.create_dataset("n_photo_cs",data = len(self.cs_data))

        for cs in self.cs_data:
            cbase = "cs_" + str(cs.id) + "/"
            cs_group.create_group(cbase)
            cs_group.create_dataset(cbase+"type",data = cs.type)
            cs_group.create_dataset(cbase+"configuration",data = cs.configuration)
            cs_group.create_dataset(cbase+"E",data = cs.E)
            cs_group.create_dataset(cbase+"sigma",data = cs.s)
            cs_group.create_dataset(cbase+"n_pts",data = len(cs.E),dtype='i')

        # write collisional data
        col_group = ion_group.create_group("collisional_data")
        cnt = 0
        line_col_id = np.full(self.nlines,-1)
        for col in self.col_data:

            # find line ID
            line_id = -1
            for i in range(self.nlines):
                if (self.line_l[i] == col.ll and self.line_u[i] == col.lu):
                    line_id = i
                    line_col_id[i] = cnt
                    break
            if (line_id < 0):
                print('no line associated with transition ', col.ll+1, col.lu+1)
                continue

            cbase = "col_" + str(cnt) + "/"
            col_group.create_group(cbase)
            col_group.create_dataset(cbase+"type",data = col.type)
            col_group.create_dataset(cbase+"Omega",data = col.O)
            col_group.create_dataset(cbase+"lev_l",data = col.ll)
            col_group.create_dataset(cbase+"lev_u",data = col.lu)
            col_group.create_dataset(cbase+"line_id",data=line_id)
            col_group.create_dataset(cbase+"T",data = col.T)
            col_group.create_dataset(cbase+"n_pts",data = len(col.T),dtype='i')

            cnt += 1

        # write index of collision id
        ion_group.create_dataset("line_col_id", data=line_col_id)

        h5f.close()


def write_simple_ion(fname,species,ion,chi,E,g):

    h5f = h5py.File(fname, 'a')

    if (not str(species) in h5f):
        species_group = h5f.create_group(str(species))
    ion_group = h5f.create_group(str(species) +'/' + str(ion))

    ion_group.create_dataset('data_source',data = np.string_("by hand"))

    ion_group.attrs['n_levels'] = len(E)
    ion_group.attrs['n_lines']  = 0

    ion_group.create_dataset('n_levels',data=int(len(E)))
    ion_group.create_dataset('n_lines',data=int(0))

            #ion_group.cattrs['ion_chi']  = self.ion_chi
    ion_group.create_dataset('ion_chi',data=chi)
    #    ion_group.create_dataset("level_i", data=self.lev)
    ion_group.create_dataset("level_g", data=g)
    ion_group.create_dataset("level_E", data=E)
    ion_group.create_dataset("level_cs", data=[0],dtype='i')


    # create fake photoionization cross-section data
    cs = CrossSection(0)
    cs.id = 0
    cs.type = "ME"
    cs.configuration = "ground"
    cs.E = np.arange(chi,100*chi,chi/10.0)
    cs.s = 6e-18*(cs.E/chi)**(-2.5)

    # write photoionization data
    cs_group = ion_group.create_group("photoion_data")
    cs_group.create_dataset("n_photo_cs",data = 1)
    cbase = "cs_" + str(cs.id) + "/"
    cs_group.create_group(cbase)
    cs_group.create_dataset(cbase+"type",data = cs.type)
    cs_group.create_dataset(cbase+"configuration",data = cs.configuration)
    cs_group.create_dataset(cbase+"E",data = cs.E)
    cs_group.create_dataset(cbase+"sigma",data = cs.s)
    cs_group.create_dataset(cbase+"n_pts",data = len(cs.E),dtype='i')


    h5f.close()


def write_nonthermal_cs(fname,datafile):

    h5f = h5py.File(fname, 'a')
    nt_group = h5f.create_group("non_thermal_cs")

    nt_group.create_dataset("source",data ="AR92")
    Z,E,chi,A,B,C,D = np.loadtxt(datafile,unpack=True)
    nt_group.create_dataset("n_data",data =len(Z))
    nt_group.create_dataset("Z",data = Z)
    nt_group.create_dataset("ion",data = Z-E)
    nt_group.create_dataset("chi",data = chi)
    nt_group.create_dataset("A",data = A)
    nt_group.create_dataset("B",data = B)
    nt_group.create_dataset("C",data = C)
    nt_group.create_dataset("D",data = D)
    h5f.close()



if __name__ == '__main__':




    # fake up simple hydrogen for tests
#    outname = 'fake_H.hdf5'
#    chi = 13.6
#    E = [0.0]
#    g = [2.0]
#    write_simple_ion(outname,1,0,chi,E,g)

    # write version
#    h5f = h5py.File(outname, 'a')
#    h5f.create_dataset("file_version",data = 1,dtype='i')
#    h5f.close()
#    exit(0)

    outname = 'cmfgen_newdata.hdf5'

#    Can redefine data sources here if we want to just
#    make a file with just few species
#    data_sources = {
#         '1':
#               {1: 'HYD/I/5dec96/hi_osc.dat'},
#         '2':
#               {1: 'HE/I/15jul15/hei_osc',
#                2: 'HE/II/5dec96/he2_osc.dat'}
#    }


    # loop over species
    for species, species_data in iteritems(data_sources):
        n_ions = 0
        # loop over ions
        for ion, ion_data_file in iteritems(species_data):
#            if (species != '8'): continue
#            if (ion-1 != 0): continue
            print(species,ion-1,ion_data_file)

            s = CMFGENDataReader(base_dir,ion_data_file,species,ion-1)
            s.write_to_file(outname,species,ion-1)
            n_ions += 1


    # cmfgen is missing data for CoI -- add in stub
    chi = 7.88101
    E   = [0.0]
    g   = [10.0]
    write_simple_ion(outname,27,0,chi,E,g)

    # cmfgen is missing data for NiI -- add in stub
    chi = 7.639878
    E   = [0.0]
    g   = [10.0]
    write_simple_ion(outname,28,0,chi,E,g)

    write_nonthermal_cs(outname,"sigma_ion_AR92.dat")

    # write version
    h5f = h5py.File(outname, 'a')
    h5f.create_dataset("file_version",data = 1,dtype='i')
    h5f.close()
