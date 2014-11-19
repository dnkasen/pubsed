#!/usr/bin/env python
import os
from matplotlib.backends.backend_pdf import PdfPages
import importlib
import time


########################################
# runs a series of tests of sedona
# output figures are put into
# test_results_%date.pdf
########################################

## executable directory and file
exec_dir   = "../src/EXEC/"
executable = "gomc"

## plot to screen
plotup = False

## sedona main codes tests to run
core_emit_into_vacuum = False
lucy_test             = False

## sedona modules tests to run
test_opacity          = True

########################################

print "TESTING EXECUTABLE: " + executable

outfile = 'test_results_' + time.strftime("%m-%d-%y") + '.pdf'
pdf = PdfPages(outfile)
print 'OUTPUT WRITTEN TO: ' + outfile

#############################################
# Script to run a test in directory: direc
#############################################
def run_one(direc,efile):


    print "------------------------------------"
    print "----- RUNNING: " + direc + "------"
    print "------------------------------------"
    print "\n\n"
    os.system("cp " + exec_dir + efile + " " + direc)
    os.chdir(direc)
    runtest = importlib.import_module(direc + ".run_test") 
    runtest.run(pdf,plotup)
    os.chdir("../")
#############################################


if (test_opacity): run_one("test_opacity","test_opacity")
if (core_emit_into_vacuum): run_one("core_emit_into_vacuum",executable)
if (lucy_test): run_one("lucy_sn",executable)

pdf.close()

