import os
import compare

def run(pdf,plotup,runcommand):

    os.system(runcommand)
    status = compare.compare(pdf)
    return status
       
