"""
AUTHOR: Chelsea E Harris
DATE: 10/03/18

PURPOSE:
This script is meant as an example for developing a new model re-formatting 
script. These scripts should be able to be called from command line OR 
be used as part of another script (via import). 
"""
import numpy as np

def read_original(rfn, read_kwargs={}):
    """
    A function to read the original model file

    INPUTS
    rfn      : [str]  path to the file
    read_args: [dict] optional arguments for reading in different ways
    """
    # This is a function that imports (needed) data from original file.
    # It may be that there are multiple ways you want to implement this, 
    # which is why there's a 
    return 0


def conversion(orig, conversion_kwargs):
    """
    A function to structure the proper array 

    INPUTS
    rfn      : [str]  path to the file
    read_args: [dict] optional arguments for reading in different ways
    """
    return 0


def write_new(wfn, write_kwargs={}):
    """
    A function to write the new model file. 
    Separated out in case you want to do checks for overwriting, 
    or make a header for using np.savetxt(), or things like that.

    INPUTS
    rfn      : [str]  path to the file
    read_args: [dict] optional arguments for reading in different ways
    """
    return 0


def main(args):
    orig = read_original(args.rfn)
    new = conversion(orig)
    write_new(new)


# This 'if' statement executes when this script is run from command line, 
# i.e., from command line,
#     >  python dev_example.py
# but does _not_ execute if the script is imported, 
# i.e., within another script,
#     from dev_example import conversion1
if __name__=='__main__':
    import argparse as ap

    parser = ap.ArgumentParser(description="Specify how to convert the code's output into Sedona input.")
    # In general, an input and output file name will be needed
    parser.add_argument('rfn', metavar='readpath', type=str, 
                        help='original model (code output) file path')
    parser.add_argument('--wfn', metavar='writepath', type=str, default='./dev_example_out.mod', 
                        help='desired new model (sedona input) file path')
    # Define the command-line arguments specific to this code:
    #parser.add_argument()

    args = parser.parse_args()

    main(args)
