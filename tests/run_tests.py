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

## executable
executable = "../src/EXEC/gomc"

## plot to screen
plotup = False

## tests to run
core_emit_into_vacuum = True
lucy_test             = True

########################################

print "TESTING EXECUTABLE: " + executable

outfile = 'test_results_' + time.strftime("%m-%d-%y") + '.pdf'
pdf = PdfPages(outfile)
print 'OUTPUT WRITTEN TO: ' + outfile

#############################################
# Script to run a test in directory: direc
#############################################
def run_one(direc):

    os.system("cp " + executable + " " + direc)
    os.chdir(direc)
    runtest = importlib.import_module(direc + ".run_test") 
    runtest.run(pdf,plotup)
    os.chdir("../")
#############################################

# --------------------------------------
# Core into vacuum test
# -------------------------------------
if (core_emit_into_vacuum):

    print "------------------------------------"
    print "----- RUNNING: emit core test ------"
    print "----- takes about 2 minutes --------"
    print "------------------------------------"
    print "\n\n"
    run_one("core_emit_into_vacuum")
########################################


# -----------------------------------------
# Lucy supernova test
# ----------------------------------------
if (lucy_test):

    print "------------------------------------"
    print "----- RUNNING: Lucy test ------"
    print "----- takes about 2 minutes --------"
    print "------------------------------------"
    print "\n\n"
    run_one("lucy_sn")
########################################


pdf.close()

