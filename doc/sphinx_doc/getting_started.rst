
=================
Getting Started
=================

-------------------
Installing the Code
-------------------

* Set the environment variable  ``SEDONA_HOME`` to point to the base directory of sedona. In bash, for example, you can add to your ``.bash_profile`` or ``.bashrc`` the line::

	export SEDONA_HOME=/Users/kasen/sedona6/


* Install lua. Source code can be found in the ``src/external/`` directory

* Download and install the Gnu Scientific Library (gsl) available at https://www.gnu.org/software/gsl/

* Download and install hdf5 from https://support.hdfgroup.org/downloads/index.html

* Set the environment variables ``GSL_DIR``, ``LUA_DIR``, and ``HDF5_DIR``. ``GSL_DIR`` and ``HDF5_DIR`` should be set automatically if you're on a supercomputer and load GSL and HDF5 modules. ``LUA_DIR`` needs to be set manually and should point to the root lua directory. For example, for me it's set to ``$SEDONA_HOME/src/external/lua-5.2.1``.
* Makefiles are in the ``makefiles`` directory. If your machine doesn't have a ``Makefile.machine`` file, use ``Makefile.general`` as a template. You'll need to set ``CXX`` to the C++ compiler you want to use, and ``CXXFLAGS`` to the compiler flags.
* Run ``chmod +x ./install.sh`` to make ``install.sh`` into something runnable.
* Compile the code by running the command ``./install.sh MACHINE`` where ``MACHINE`` corresponds to the ``Makefile.MACHINE`` you want to use. This will copy all source files into ``build`` and compile.
* Other options for ``./install.sh`` are:
  * ``help``: prints usage statement
  * ``clean``: deletes source, object, and binary files in ``build``
  * ``realclean``: deletes entire ``build`` directory

* If compilation is successful, the executable file ``sedona6.ex`` will appear in the ``src/`` directory. Copy this to the directory where you would like to run the code.



Python is currently used for plotting and testing scripts, but is not needed to build and run |sedona| itself (NB: Python may
eventually replace lua for parameter files). Python 2.7 is being used. On linux to get the necessary packages try::

	sudo apt-get install python-numpy python-scipy python-matplotlib python h5py



-------------------
Running the Code
-------------------

To get started with |sedona|, you can try running a simple transport calculation. Copy the executable file |sedona.ex| into a directory where you want to run, for example type::

  cp sedona6.ex examples/simple_examples/spherical_lightbulb/1D/

Change to that directory, and run the code by typing::

	./sedona6.ex param.lua


where *param.lua* is the file name of the run time parameter file (if no file name is given, the name *param.lua* is assumed).  See :ref:`parameter_files` for description of setting runtime parameters. At minimum, the parameter file must point to three files


* A model file  (an ascii .mod file for 1D runs or a hdf5 file for multi-D) giving the physical conditions (e.g., density, composition) of the setup. See :ref:`model_files` for description of formats.

* An atomic data file in hdf5 format. Such files are available in the ``data`` directory, with data compiled from various sources.

* A  defaults file, giving the default settings for all runtime parameters. The standard is ``defaults/sedona\_defaults.lua``, although users can point to their own modified defaults file.


The code generates several files. The *plt_.?????.h5* files contain data in hdf5 format describing the grid properties (e.g., density, temperature, radiation field). If hdf5 tools are installed, type::

	h5ls plt_00001.h5

to see the file contents. For 1D models, a *plt_?????.dat* gives some of this data in ascii format.

The code also generates a  *spectrum_?.h5* hdf5 file (and for 1D models a *spectrum_?.dat* asci file) giving the output spectrum.

You will find in the ``tools/`` directory a python library file *sedona.py* that provides functions to read and analyze the data in the plot and spectrum files (NB: this needs to be developed)
