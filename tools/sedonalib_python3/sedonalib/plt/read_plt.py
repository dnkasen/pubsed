from .SedonaPltFile import *

def read_pltfile(name):

    if (isinstance(name, int)):
        name = str(name)
        name = "plt_" + name.rjust(5, '0') + ".h5"

    return SedonaPltFile(name)