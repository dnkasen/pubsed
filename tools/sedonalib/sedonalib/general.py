from spectrum import *

def read(fname):

     ascii_extensions = [".dat",".mod",".txt"]
     hdf5_extensions = [".h5",".hdf5"]

     fileformat = None
     filetype   = None

     # create models from an ascii file format
     if any(fname.endswith(x) for x in ascii_extensions):
         fileformat = "ascii"

     # create models from an hdf5 file format
     if any(fname.endswith(x) for x in hdf5_extensions):
         fileformat = "hdf5"

         import h5py
         fin = h5py.File(fname,'r')

         if ("nu/" in fin):
             return read_spectrum(fname)


     # throw error if file format unknown
     if (fileformat is None):
         message = "Unknown file type, allowed values are "
         message += str(ascii_extensions + hdf5_extensions)
         raise Exception(message)
