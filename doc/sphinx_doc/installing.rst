=================
Installing Sedona
=================

-------------------
Getting the Code
-------------------

The code is available for download on github `here <https://github.com/dnkasen/sedona6>`_.
Or you can clone from the command line using::

  git clone git@github.com:dnkasen/sedona6.git

.. _dependencies:

--------------------------
Installing Dependencies
--------------------------

[sedona] requires a C++ compiler for instalation. In addition. the following dependencies must be installed:

* `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ (file formating)
* `GSL <https://www.gnu.org/software/gsl/doc/html/>`_ (gnu scientific library)
* `lua <http://www.lua.org>`_ (scripting language used for parameter files)

[sedona] is parallelized using a hybrid of MPI and openMP. To run in parallel, you must also have installed an MPI  and/or an openMP distribution.

These packages may already been installed on clusters and computing systems. If not, they can be conveniently installed using a package manager -- see below for details.

If your dependencies are not installed in a standard location (e.g., usr/local/) then you will need to set the ``GSL_DIR, HDF5_DIR``, and ``LUA_DIR`` environment variables to specify the path to the instalations. In the bash shell, for example, add to your bash.profile::

  export  GSL_DIR=/base/path/of/gsl/
  export HDF5_DIR=/base/path/of/hdf5/
  export  LUA_DIR=/base/path/of/lua/

where the ``/base/path/of/xxx`` should be replaced with the full base path where the package is installed.


A python distribution is not required to run [sedona] itself, however it can be useful for  working with input and output files (in particular the `h5py <https://www.h5py.org>`_ package allows for easy creation and reading of hdf5 files).  The [sedona] distribution comes with a python package called *sedonalib* that contains useful tools for
generating, reading, and plotting files. The `anaconda <https://www.anaconda.com/products/individual>`_ python distribution  makes it relatively easy to install the needed python packages.


....

**Installing Dependencies on a Mac**



On a mac, the *homebrew* package manager provides a convenient way to install dependences. Type ``brew help`` to check if homebrew it is already installed on your machine. If the command is not found, it can be installed by issuing
the command given on the
`homebrew installation page <https://docs.brew.sh/Installation>`_

Once installed, the necessary [sedona] dependencies can be
installed with the following commands::

  brew install lua
  brew install gsl
  brew install hdf5

If you wish to run mpi parallel jobs, you should also install mpich::

  brew install mpich

If you lack a C++ compiler, you can install the gcc compiler through ``brew install gcc``. Alternatively, mac Mac provides a C++ compiler through its Xcode tools.


....

**Installing Dependencies on a Linux System**

On Debian, Ubunto and some other linux systems, the *apt* package manager provides a convenient way to install dependences. If *apt* is installed on your machine, you can update it using::

  sudo apt-get update

Root privileges are required. You can then install the necessary dependencies using::

  sudo apt-get install gsl-bin
  sudo apt-get -s install lua5.2
  sudo apt-get install libhdf5-dev

If you wish to run mpi and openMP parallel jobs, you should also install::

  sudo apt-get install  mpich
  sudo apt install libomp-dev

....

**Installing Dependencies on a Cluster**

Many clusters and supercomputer centers will have some of the dependencies already installed.  Use the module load





-----------------------
Compiling the Code
-----------------------


To compile [sedona], first change to the ``src/`` directory in the [sedona] distribution, and type::

  chmod +x ./install.sh

to ensure that the compile script is executable.  [sedona] comes with a set of makefiles to facilitate compilation on different machines.  To see a list of machines that have makefiles, type::

  ./install.sh help

If you see a machine name appropriate to your system, you can attempt compiling the code by typing::

  ./install.sh <MACHINE>

where ``<MACHINE>`` is the name of the machine.

If you don't see a relevant Makefile for your machine, or if compilation fails, you should modify a Makefile directly. Open the ``makefiles/Makefile.general`` file and make sure that the ``CXX`` variable is set to your  C++ compiler. You can also change the ``CXX_FLAGS`` variables to set the compiler flags of your choosing. You can also check that all dependent packages are installed and that the GSL_DIR, HDF5_DIR and LUA_DIR environment variables are set,  as described in :ref:`dependencies`


Once compilation is successful, the executable file ``sedona6.ex`` will appear in the ``src/`` directory. Copy this to the directory where you would like to run the code.



-------------------------
Installing Python Tools
--------------------------
