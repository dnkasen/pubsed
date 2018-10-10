"""
Author: Chelsea E Harris
Creation Date: Oct 1, 2018
Purpose: Convert Kepler output into a Sedona input model for light curve and spectra creation.
"""

import numpy as np

# Maps between the way Kepler refers to isotopes (e.g., 'h1') and the way Sedona does (e.g., '1.1')
Sed2Kep = {'1.1':'h1', '2.4': 'he4', '6.12':'c12', '7.14': 'n14', '8.16':'o16', '10.20':'ne20',
           '12.24':'mg24', '14.28':'si28', '16.32':'s32', '18.36':'ar36', '20.40':'ca40', 
           '22.44':'ti44', '24.48':'cr48', '26.52':'fe52', '26.54':'fe54', '26.56':'fe56', 
           '27.56':'co56', '28.56':'ni56'}
Kep2Sed = {}
for key in Sed2Kep.keys():
    Kep2Sed[Sed2Kep[key]] = key
# These species are in Kepler but not in Sedona (or for some other reason 
# it is desired to lump them in with another element).
# In the conversion process they will be lumped with other elements
Kep2Sed['nt1'] = '' # to be lumped with 1.1
Kep2Sed['pn1'] = '' # ... 1.1
Kep2Sed['he3'] = '' # ... 2.4
Kep2Sed['fe54'] = '' # ... 26.52

def read_original(rfn, verbose=False):
    """
    Reads a Kepler table into a dictionary. Kepler output tables are 
    a difficult format to read in because there aren't an equal number of
    columns per row, nor a specified header, and there can be multiple tables
    on one row (in the "header")
    INPUTS
    rfn : path to the Kepler model (".link") file
    """
    if verbose: print("Reading file.")

    # Because the file is all weird, I'm just opening it the non-numpy way
    with open(rfn,'r') as rf:
        all_lines = rf.read().split('\n')

    # Separate the top header from the rest, and locate the "blocks" of data
    mat = []
    i_block_start = []
    n = 0
    main_header = ''
    current_cols = ''
    for line in all_lines:
        cols = line.split()

        # Skip empty or "1:" lines
        if len(cols) < 2: 
            continue

        if cols[0] == 'zone': 
            # sometimes the blocks are divided up and they get a new
            # column labeling line each time
            if line != current_cols:
                i_block_start.append(n)
                current_cols = line
            else:
                continue
        # skip summation lines
        elif cols[0] == 'total':
            continue
        # stop if we've reached the end of the useful info
        elif cols[0] == 'parameter':
            break

        # add columns to the matrix -- keep the hader in, so that
        # the line number of the block starting can be used 
        mat.append(cols)
        # if we haven't started the blocks yet, we're still in the header
        if len(i_block_start)==0: 
            main_header += line

        n+=1

    # Add fake block start so we can do the next for-loops
    i_block_start.append(len(all_lines))
    # Put blocks into arrays, taking out the "zone" column
    blocks = [ np.array(mat[i_block_start[0]+1: i_block_start[0+1]])[:,1:] ]
    headers = [ mat[i_block_start[0]][1:] ]
    # also take out mass column for blocks after the first one
    blocks += [ np.array(mat[i_block_start[i]+1: i_block_start[i+1]])[:,2:]  for i in range(1,len(i_block_start)-1) ]
    headers += [ mat[i][2:] for i in i_block_start[1:-1] ]

    # Truncate the second block because its 'burn off' entries are counted as two columns! 
    # And I don't need them
    blocks [1] = blocks [1][:,:9]
    headers[1] = headers[1][:9]

    # Check that the blocks are OK in shape
    print(len(blocks))
    for i_block, block in enumerate(blocks[1:]):
        if block.shape[0] - blocks[0].shape[0] != 0: 
            print('FATAL ERROR: It seems that the script is not identifying blocks correctly.')
            print('  ... got size {} instead of {} in Block {}'.format(block.shape[0], blocks[0].shape[0], i_block+2))
            return -1

    # Merge the blocks into one table, same with the column headers

    # populate the headers array
    all_hdrs = []
    species = []
    for i, hdr_line in enumerate(headers):
        all_hdrs += hdr_line
        if i>=2:
            species += hdr_line
    all_hdrs = np.array(all_hdrs)
    species = species[:-1]

    # the data array
    n_rows = blocks[0].shape[0]
    n_cols = sum( [block.shape[1] for block in blocks] )
    
    tab = np.zeros((n_rows, n_cols)) # each block is contributing a "zone" column
    j0 = 0
    for block in blocks:
        try:
            tab[:, j0:j0+block.shape[1]] = block.astype(float)
        except ValueError:
            for j in range(block.shape[1]):
                try:
                    tab[:, j0+j] = block[:,j].astype(float)
                except ValueError:
                    if verbose:
                        print('Setting the {} column to zero because it isn\'t a float'.format(all_hdrs[j0+j]))
                        print('  e.g., value of '+repr(block[0,j]))
        j0 += block.shape[1]

    # Make into a dictionary
    Table = {}
    for j in range(all_hdrs.size):
        Table[all_hdrs[j]] = tab[:,j]

    return Table, main_header, species


def read_stripped(rfn, time, verbose=False):
    """
    Reads a stripped version (no strings) of a  Kepler table into a dictionary.
    INPUTS
    rfn : path to the Kepler model file
    time : simulation time (seconds)
    """
    blk1 = ['mass', 'vmass', 'radius', 'density', 'temp', 'energy', 'sdot', 'snuc', 'sburn','abar', 'lumin', 'velocity']
    blk2 = ['entropy', 'mixtime', 'eta', 'pressure', 'stress', 'opacity', 'ye', 'yeburn', 'limnuc(cyczb)', 'dtburn']
    blk3 = [ 'nt1', 'h1', 'pn1', 'he3', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24'] 
    blk4 = ['si28', 's32', 'ar36', 'ca40', 'ti44', 'cr48', 'fe52', 'fe54', 'ni56', '1-tot']
               

    zones = np.loadtxt(rfn, usecols=[0], dtype=int)
    N_zones = np.unique(zones).size

    with open(rfn,'r') as rf:
        all_lines = rf.read().split('\n')

    Table = {}

    Table['vmass'   ] = np.array( [ float(line.split()[ 2])    for line in all_lines[:N_zones] ] )
    Table['radius'  ] = np.array( [ float(line.split()[ 3])    for line in all_lines[:N_zones] ] )
    Table['density' ] = np.array( [ float(line.split()[ 4])    for line in all_lines[:N_zones] ] )
    Table['temp'    ] = np.array( [ float(line.split()[ 5])    for line in all_lines[:N_zones] ] )
    Table['velocity'] = np.array( [ float(line.split()[-1])    for line in all_lines[:N_zones] ] )

    species = blk3 + blk4[:-1]

    x  = np.zeros((N_zones,19))
    for i in range(N_zones):
        x[i,:10] = [ float(col) for col in all_lines[  N_zones+i].split()[2:  ] ]
        x[i,10:] = [ float(col) for col in all_lines[2*N_zones+i].split()[2:-1] ]

    for j, spec in enumerate(species):
        Table[spec] = x[:,j]
    
    return Table, '', species


def convert(Table, header, species, 
            time=0, verbose=False):
    """
    This function takes the information from the Kepler table and generates
    the table that Sedona needs.
    It does not generate the header that Sedona needs, but does read (and return) 
    the time.
    Note that the 'nt1' and 'pn1' columns are lumped into '1.1' with 'nt1' multiplied
    by 0.1,
    the 54-Co mass fraction is set to 1e-10, 
    and 52-Fe includes 'fe52' and 'fe54'    

    INPUTS
    Table  : [dict]  key, value pairs are column names and columns from the Kepler blocks
                     (output of read_original())
    header : [str]   the Kepler model file before the blocks start
                     (output of read_original())
    species: [str[]] list of the species to put into the Sedona model, written as Sedona
                     types (e.g. '1.1')
                     (output of read_original() or specified)
    verbose: [bool]  whether to print messages

    OUTPUTS
    mat    : [np.ndarray(float)] the main body of the Sedona model
    time   : [float]             model time (sec)
    """

    if verbose:
        print('Converting data to Sedona version')

    time = float(header.split()[10]) if time==0 else time
    n_zones = len(Table[list(Table.keys())[0]])

    # number of gas properties needed by Sedona
    N_needed = 3
    # initialize the matrix that will be ouptut as the model
    mat = np.zeros((n_zones, N_needed+len(species)))

    mat[:,0] = Table['radius']/time
    mat[:,1] = Table['density']
    mat[:,2] = Table['temp']

    for i, s in enumerate(species):
        # If the species list is coming from Kepler, some 
        # elements will have been unknown
        if s=='':
            continue
        elif s=='1.1':
            x = Table['h1'] + Table['nt1'] + Table['pn1']
        elif s=='2.4':
            x = Table['he3'] + Table['he4']
        elif s=='27.54' or s=='27.56':
            x = 1e-10*np.ones(n_zones)
            print('Setting Co-54 to 1e-10')
        elif s=='26.52' or s=='26.56':
            x = Table['fe52'] + Table['fe54']
        else:
            try:
                x = Table[Sed2Kep[s]]
            except KeyError:
                print(s, Se)

        mat[:, i+N_needed] = x

    return mat, time


def overwrite_safety(wfn):
    """
    This function checks if wfn already exists, and re-names the
    new file destination if so.
    """
    import glob
    
    exists = len(glob.glob(wfn))
    if exists:
        print('Will not write file '+wfn+' -- file exists!')
        base = wfn
        i=1
        while exists:
            wfn = base+'.'+str(i)
            exists = len(glob.glob(wfn))
            i+=1
    print('Writing to '+wfn)

    return wfn


def write_new(wfn, time, species, sim_data, rfn='', 
              grid_type='1D_sphere', hydro_type='SNR', overwrite=False):
    """
    A function to generate the Sedona header and write the sedona model file

    INPUTS
    wfn  : [str] path for new file
    time : [float] time of model in seconds
    species: [list(str)] species to include, in Sedona notation (e.g., '1.1')
    sim_data: [np.ndarray(float)] main Sedona body (output of convert())
    """
    # Check that we have a write-file name/path:
    if wfn=='':
        if rfn=='':
            wfn = './KeplerOut.mod'
        else:
            wfn = rfn.replace(rfn.split('.')[-1], 'mod')
    # If safety is desired, check this name
    if not overwrite:
        wfn = overwrite_safety(wfn)

    # Write the Sedona header
    hdr = grid_type + '  ' + hydro_type + '\n'
    hdr += '{:d} 0.0 {:.5e} {:d}\n'.format(sim_data.shape[0], time, len(species))
    hdr += ' '.join(species)

    # store it all!
    np.savetxt(wfn, sim_data, header=hdr, comments='')
    
            

def main(args):
    # Read in the Kepler file
    if not args.stripped_file:
        readout = read_original(args.rfn, verbose=args.verbose)
    else:
        readout = read_stripped(args.rfn, args.time, verbose=args.verbose)
    if type(readout)==int: 
        return 0

    Table, header, species = readout
    # Turn it into a data matrix for Sedona
    spec2use = args.specs if len(args.specs)>0 else [ Kep2Sed[s] for s in species ]
    mat, time = convert(Table, header, spec2use, time=args.time, verbose=args.verbose)
    # Write to new file
    write_new(args.wfn, time, spec2use, mat, 
              rfn=args.rfn, grid_type=args.grid, hydro_type=args.hydro, overwrite=args.overwrite)



if __name__=='__main__':
    import argparse

    parser = argparse.ArgumentParser()
    # Make arguments for this converter
    parser.add_argument('rfn', type=str, help='path to file')

    parser.add_argument('--time', type=float, default=0.0, help='time of the simulation in seconds')

    parser.add_argument('-sf', '--stripped_file', action='store_true', 
                        help='whether the file represents raw Kepler output or stripped-down version')

    parser.add_argument('--wfn', type=str, default='', help='path for desired output file')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='(flag) print messages to screen during execution')
    parser.add_argument('--specs', type=str, nargs='+', default=[], 
                        help='list of species for sedona (default is all that are in the Kepler file) formatted like --specs = s1 s2 s3 ...')

    parser.add_argument('--grid', type=str, default='1D_sphere',
                        help='sedona grid parameter to use')

    parser.add_argument('--hydro', type=str, default='SNR',
                        help='sedona hydro parameter to use')

    parser.add_argument('--overwrite', action='store_true', 
                        help='(flag) throw caution to the wind when writing the new file')

    # Use the arguments to make the appropriate conversion
    args = parser.parse_args()
    main(args)
