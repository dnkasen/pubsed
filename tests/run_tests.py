#!/usr/bin/env python
import os

## executable
executable = "../src/exec/gomc"
print "TESTING EXECUTABLE: " + executable


########################################
print "------------------------------------"
print "----- RUNNING: emit core test ------"
print "----- takes about 2 minutes --------"
print "------------------------------------"
print "\n\n"

dir = "core_emit_into_vacuum"
os.system("cp " + executable + " " + dir)
os.chdir(dir)
os.system("./gomc")
os.system("python compare.py")
os.system("rm out_optical_1.spec ray_* gomc")
os.chdir("../")
########################################
