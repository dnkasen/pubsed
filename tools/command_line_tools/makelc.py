#!/usr/bin/env python
#
#
#
# Instructions:
# --------------
#   To run, use the command "python lcfilt.py -s <input file> -b <band1,band2,...>
#   where <input file> is the path to the spectrum.h5 file
#   and <band1,band2,...> is a list of comma-separated filters
#
#   Example Usage: python lcfilt.py -s /Users/dkhatami/data/spectrum.h5 -b U,B,V
#       ...generates U,B,V AB magnitude light curves from the spectrum file
#
#   For a list of filters, use "python lcfilt.py --bands"
# -------------
#   Return: File called "lightcurve.out" in the format of
#   	Columns: 1. Time (Days) 2. Bolometric Luminosity (erg/s) 3. Bol. Mag. [4,5,...] Filter magnitudes
# 	 	Invalid elements are assigned a value of "0"
#		Floors the magnitudes to 0 if <0
# -------------
#
# Provided filter data taken from Charlie Conroy's Flexible Stellar Population Synthesis (FSPS) Code
#  See the doc/manual.pdf in the FSPS github for how the magnitudes are defined
#

import sys, getopt
import optparse
import h5py
import numpy as np
import sedonalib as sed

def main(argv):

    parser = optparse.OptionParser(usage="%prog <spectrum_filename> -b <bands> [options]")
    parser.add_option("-f","--filters",action="store_true", dest="show_filters", help="Print out list of all available filter bands")
    parser.add_option("-b","--bands",dest="bands",help="List of bands to calculate, default=bolometric")
    parser.add_option("-p","--plot",action="store_true",dest="plotit",help="Plot up light curves to screen")
    parser.add_option("-o","--outfile",dest="outfile",help="Name of outputfile, default=\"lightcurve.out\"")

    (opts, args) = parser.parse_args()

    if (opts.show_filters):
        filt = sed.Filter()
        bands = tuple(filt.Bands)
        print('-----------------------------')
        print('List of Available Filters:')
        print('-----------------------------')
        for b in bands:
            print(b)
        print('-----------------------------')
        print('See file FILTER_LIST for details/references')
        print('-----------------------------')
        sys.exit(0)

    if (len(args) < 1):
        print('ERROR: Need to supply spectrum filename\n')
        parser.print_help()
 #usage is: \n  makelc.py <spectrum.h5> -b <band1,band2,...>')
        exit(1)
    else:
        specfile = args[0]

    if (opts.bands):
        bandlist = opts.bands.split(',')
    else:
        bandlist = ['B','V','Cousins_R']

    if (opts.outfile):
        outfile = opts.outfile
    else:
        outfile = "lightcurve.out"

    print("Calculating light curves for bands: ",bandlist)
    spec = sed.read_spectrum_file(specfile)
    spec.write_lightcurve(band=bandlist,magnitudes='AB',outfile=outfile)

    if (opts.plotit):
        spec.plot_lightcurve(band=bandlist,magnitudes='AB')

    print("Output written to {}".format(outfile))

if __name__ == "__main__":
    main(sys.argv[1:])
