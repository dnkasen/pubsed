"""
Author: Chelsea E Harris
Creation Date: Oct 1, 2018
Purpose: Convert Kepler output into a Sedona input model for light curve and spectra creation.
"""

import argparse # Tool for parsing command line arguments
import numpy

def read_original(rfn):
    


def main(args):



if __name__=='__main__':
    parser = argpase.ArgumentParser()
    # Make arguments for this converter
    parser.add_argument('rfn', type=str, help='path to file')

    parser.add_argument('time', type=float, help='time of the simulation in seconds')

    parser.add_argument('--wfn', type=str, help='path for desired output file')

    # Use the arguments to make the appropriate conversion
    args = parser.parse_args()
    main(args)
