from SedonaModel import *
from Sedona1DSphereModel import *

def read_model(name,type=type):

    ascii_extensions = [".dat",".mod",".txt"]
    hdf5_extensions = [".h5",".hdf5"]

    fileformat = ""

    # create models from an ascii file format
    if any(name.endswith(x) for x in ascii_extensions):
        fileformat = "ascii"
        with open(name) as fin:
            geometry,type = fin.readline().split()

            if (geometry == "1D_sphere"):
                return Sedona1DSphereModel(name)
            if (geometry == "2D_cyln"):
                return Sedona2DCylnModel(name)

    # create models from an hdf5 file format
    if any(name.endswith(x) for x in hdf5_extensions):
        fileformat = "hdf5"

    # through error if file format unknown
    if (fileformat == ""):
        message = "Unknown filename extension, allowed values are "
        message += str(ascii_extensions + hdf5_extensions)
        raise Exception(message)




def new_model(geometry,dims=None,type=type,rout=None):

    if (geometry == "1D_sphere"):
        return Sedona1DSphereModel(dims=dims,rout=rout)
