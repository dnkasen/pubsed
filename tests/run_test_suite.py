#!/usr/bin/env python
import os, sys
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import time
import optparse

########################################
# runs a series of tests of sedona
# output figures are put into
# test_results_%date.pdf
########################################

parser = optparse.OptionParser()
parser.add_option("-n",dest="nproc")
(opts, args) = parser.parse_args()

nproc = 1
if (opts.nproc): nproc = opts.nproc

## executable directory and file
exec_dir   = "../src/EXEC/"
executable = "gomc"
runcommand = "mpirun -np " + str(nproc) + " ./gomc"

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

    homedir = os.getcwd()

    print "------------------------------------"
    print "----- RUNNING: " + direc 
    print "------------------------------------"
    print "\n\n"
    os.system("cp " + exec_dir + efile + " " + direc)
    os.chdir(direc)
    sys.path.append(os.getcwd())
    import run_test  
    run_test.run(pdf,plotup,runcommand)
    sys.path.remove(os.getcwd())
    os.chdir(homedir)
#############################################


# loop over tests in suite_test_list file and run them
fin = open("suite_test_list","r")
for line in fin:
    line = line.rstrip('\n')
    line = line.rstrip(' ')
    if (os.path.isdir(line)): 
        run_one(line,executable)
    

pdf.close()

