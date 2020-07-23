#!/usr/bin/env python
import os, sys
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import time
import optparse
from importlib import import_module
import timeit


###############################################
# runs a series of tests of sedona
# by default, the tests to be run are read
# from the file "suite_test_list"
#
# output figures showing test results are put
# into in the pdf file:
#  test_results_%date.pdf
# (where %date is the current date)
#
# Usage:
#  python run_test_suite.py [options]
# Options:
#
#   -n 2
#   (will run on 2 mpi ranks, more generally -n nranks
#   will run on nranks ranks where nranks is an integer)
#
#  --testlist "spherical_lightbulb/1D","lucy_supernova/1D"
#   (will run the tests given after the --testlist flag)
#
#  --testfile "my_tests.txt"
#   (will read the file "my_tests.txt" or whatever filename you wante
#   to get the list of files to run)
###########################################

parser = optparse.OptionParser()
parser.add_option("-n",dest="nproc")
parser.add_option("--testlist","-t",dest="testlist",type='string')
parser.add_option("--testfile","-f",dest="testfile",type="string")
parser.add_option("--outfile","-o",dest="outfile",type="string")
parser.add_option("-v",action="store_true",dest="verbose")
parser.add_option("--mpi_args",dest="mpi_args", type="string")

(opts, args) = parser.parse_args()

nproc = 1
if (opts.nproc): nproc = opts.nproc

## executable directory and file
exec_dir   = "../src/"
executable = "sedona6.ex"
paramname = "param.lua"

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

if (opts.mpi_args):
  mpi_args = opts.mpi_args
else:
  mpi_args = ""

runcommand = "mpirun -np " + str(nproc) + " " + mpi_args + " ./" + executable + " " + paramname
if (not opts.verbose):
    runcommand = runcommand + " >> " + outfile + " 2>&1"


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
print("Testing SEDONA code on " + date + "\n")
os.system(line)

# get number of threads
try:
    nthreads = os.environ['OMP_NUM_THREADS']
except:
    nthreads = 1

print("Will run " + str(len(testlist)) + " tests on " +str(nproc)  + " mpi ranks,"),
print("with " + str(nthreads) + " threads per rank\n");
print("------------------------------------------")
i = 1
for this_test in testlist:
    print(str(i) +') ' + this_test)
    i = i + 1
print("------------------------------------------")
print("\n")


total_status = 0
cnt = 1
for this_test in testlist:

    starttime = timeit.default_timer()

    hdr = "\n\n------------------------------------\n"
    hdr = hdr + "-- RUNNING: " + this_test + "\n"
    hdr = hdr + "------------------------------------\n"
    cmd = "echo \"" + hdr + "\"" +  " >> " + outfile
    os.system(cmd)

    print("------------------------------------------")
    print("- test " + str(cnt) + ") " + this_test)

    os.system("cp " + exec_dir + executable + " " + this_test)
    os.chdir(this_test)

    # immport and run script
    sys.path.append(os.getcwd())
    import run_test
    sys.path.remove(os.getcwd())
    try:
      status = run_test.run_test(pdf,runcommand)
    except:
      status = 200
    del sys.modules['run_test']

    stoptime = timeit.default_timer()
    print("time = {:.1f} seconds".format(stoptime - starttime))

    if (status == 0):
        print("PASSED")
    elif (status == 200):
        print("CRASHED")
    else:
        print("FAILED: error code = ", status)
    print("------------------------------------------\n")
    total_status += status

    # return home
    os.chdir(homedir)
    cnt += 1


pdf.close()
