#!/usr/bin/env python

import sys, os
import optparse

parser = optparse.OptionParser()

parser.add_option("-l",dest="log")
parser.add_option("--xr",dest="xrange")
parser.add_option("--yr",dest="yrange")

(opts, args) = parser.parse_args()


#for opt,arg in opts:
#  print opt,arg

cmd = "bokeh serve --show  ./plot_lc.py --args "
for f in args:
    cmd += " " + f

os.system(cmd)


#n = len(args)
