#!/usr/bin/env python
import os, sys
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import time
import optparse
from importlib import import_module
import timeit


########################################
# runs a series of tests of sedona
# output figures are put into
# test_results_%date.pdf
########################################

parser = optparse.OptionParser()
parser.add_option("-n",dest="nproc")
parser.add_option("--testlist","-t",dest="testlist",type='string')
parser.add_option("--testfile","-f",dest="testfile",type="string")
parser.add_option("--outfile","-o",dest="outfile",type="string")
parser.add_option("-v",action="store_true",dest="verbose")

(opts, args) = parser.parse_args()

nproc = 1
if (opts.nproc): nproc = opts.nproc

## executable directory and file
exec_dir   = "../src/"
executable = "sedona6.ex"

homedir = os.getcwd()
date = time.strftime("%m-%d-%y")

# setup output files
if (opts.outfile):
    outfile = homedir + '/' + opts.outfile + '.txt'
    pdffile = homedir + '/' + opts.outfile + '.pdf'
else:
    outfile = homedir + '/test_results_' + date + '.txt'
    pdffile = homedir + '/test_results_' + date + '.pdf'
pdf = PdfPages(pdffile)


runcommand = "mpirun -np " + str(nproc) + " ./" + executable
if (not opts.verbose):
    runcommand = runcommand + " >> " + outfile


########################################

# get list of tests eihter from input flag
if (opts.testlist):
    testlist = opts.testlist.split(',')
# or else from the default list file
else:
    testfile_name = "suite_test_list"
    if (opts.testfile): testfile_name = opts.testfile
    fin = open(testfile_name,"r")
    testlist = []
    for line in fin:
        line = line.rstrip('\n')
        line = line.rstrip(' ')
        if (os.path.isdir(line)):
            testlist.append(line)

line = "echo \"Testing SEDONA code on " + date + "\"  > " + outfile
print "Testing SEDONA code on " + date + "\n"
os.system(line)
print "Will run " + str(len(testlist)) + " tests on " +str(nproc)  + " cores\n"
print "------------------------------------------"
i = 1
for this_test in testlist:
    print str(i) +') ' + this_test
    i = i + 1
print "------------------------------------------"
print "\n"


total_status = 0
cnt = 1
for this_test in testlist:

    starttime = timeit.default_timer()

    hdr = "\n\n------------------------------------\n"
    hdr = hdr + "-- RUNNING: " + this_test + "\n"
    hdr = hdr + "------------------------------------\n"
    cmd = "echo \"" + hdr + "\"" +  " >> " + outfile
    os.system(cmd)

    print "------------------------------------------"
    print "- test " + str(cnt) + ") " + this_test

    os.system("cp " + exec_dir + executable + " " + this_test)
    os.chdir(this_test)

    # immport and run script
    sys.path.append(os.getcwd())
    import run_test
    sys.path.remove(os.getcwd())
    status = run_test.run_test(pdf,runcommand)
    del sys.modules['run_test']

    stoptime = timeit.default_timer()
    print 'time = {:.1f} seconds'.format(stoptime - starttime)

    if (status == 0):
        print "PASSED"
    else:
        print "FAILED: error code = ", status
    print "------------------------------------------\n"
    total_status += status

    # return home
    os.chdir(homedir)
    cnt += 1


pdf.close()
