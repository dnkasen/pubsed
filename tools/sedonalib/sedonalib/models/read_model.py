from SedonaModel import *
from Sedona1DSphereModel import *
from Sedona2DCylnModel   import *
from Sedona3DCartModel   import *

def read_model(name,type=type):

    ascii_extensions = [".dat",".mod",".txt"]
    hdf5_extensions = [".h5",".hdf5"]

    fileformat = ""

    # create models from an ascii file format
    if any(name.endswith(x) for x in ascii_extensions):
        fileformat = "ascii"
        with open(name) as fin:
            geometry,type = fin.readline().split()

    # create models from an hdf5 file format
    if any(name.endswith(x) for x in hdf5_extensions):
        import h5py

        fileformat = "hdf5"
        fin = h5py.File(name,'r')
        ndims = len((fin['rho']).shape)
        if (ndims == 1):
            geometry = "1D_sphere"
        if (ndims == 2):
            geometry = "2D_cyln"
        if (ndims == 3):
            geometry = "3D_cart"

    # throw error if file format unknown
    if (fileformat == ""):
        message = "Unknown filename extension, allowed values are "
        message += str(ascii_extensions + hdf5_extensions)
        raise Exception(message)

    if (geometry == "1D_sphere"):
        return Sedona1DSphereModel(name)
    if (geometry == "2D_cyln"):
        return Sedona2DCylnModel(name)
    if (geometry == "3D_cart"):
        return Sedona3DCartModel(name)


def new_model(geometry,dims=None,type=type,rout=None):

    if (geometry == "1D_sphere"):
        return Sedona1DSphereModel(dims=dims,rout=rout)
    if (geometry == "2D_cyln"):
        return Sedona2DCylnModel(dims=dims,rout=rout)
    if (geometry == "3D_cart"):
        return Sedona3DCartModel(dims=dims,rout=rout)
