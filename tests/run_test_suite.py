#!/usr/bin/env python
import os, sys
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import time
import optparse
from importlib import import_module



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


########################################

# get list of tests in suite_test_list file 
fin = open("suite_test_list","r")

testlist = []
for line in fin:
    line = line.rstrip('\n')
    line = line.rstrip(' ')
    if (os.path.isdir(line)): 
        testlist.append(line)


print "TESTING EXECUTABLE: " + executable
print "Will run " + str(len(testlist)) + " tests"
print "------------------------------------------"
i = 0
for this_test in testlist:
    print str(i) +') ' + this_test
    i = i + 1
print "\n"


for this_test in testlist:
    homedir = os.getcwd()

    print "------------------------------------"
    print "----- RUNNING: " + this_test 
    print "------------------------------------"
    print "\n"

    os.system("cp " + exec_dir + executable + " " + this_test)
    os.chdir(this_test)

    # immport and run script
    sys.path.append(os.getcwd())
    import run_test
    sys.path.remove(os.getcwd())
    run_test.run_test(pdf,runcommand)
    del sys.modules['run_test']

    # return home
    os.chdir(homedir)



pdf.close()

