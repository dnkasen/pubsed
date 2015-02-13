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
runcommand = "mpirun -np 2 ./gomc"

outfile = 'test_results_' + time.strftime("%m-%d-%y") + '.pdf'
pdf = PdfPages(outfile)
print 'OUTPUT WRITTEN TO: ' + outfile
print 'run command: ' + runcommand
## plot to screen
plotup = False


########################################

print "TESTING EXECUTABLE: " + executable



#############################################
# Script to run a test in directory: direc
#############################################
def run_one(direc,efile):


    print "------------------------------------"
    print "----- RUNNING: " + direc 
    print "------------------------------------"
    print "\n\n"
    os.system("cp " + exec_dir + efile + " " + direc)
    os.chdir(direc)
    runtest = importlib.import_module(direc + ".run_test") 
    runtest.run(pdf,plotup,runcommand)
    os.chdir("../")
#############################################


fin = open("test_list","r")
for line in fin:
    line = line.rstrip('\n')
    line = line.rstrip(' ')
    if (os.path.isdir(line)):  run_one(line,executable)

pdf.close()

