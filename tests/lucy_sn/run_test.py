import os
import compare

def run(pdf,plotup,runcommand):

    ###########################################
    # run the code
    ###########################################
    os.system(runcommand)
       
    ###########################################
    # compare the output
    ###########################################
    compare.compare(pdf)

    ###########################################
    # clean results
    ###########################################
    os.system("rm optical_spectrum* gamma* ray_* gomc")
