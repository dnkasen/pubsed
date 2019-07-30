=================
Python Tools
=================

-------------------
Overview
-------------------
The python tools is a small library of tools and classes useful for handling
sedona in- and output. So far it consists of the following classes:

* ``ModelFile`` class (h5 files ``MODEL.h5`` and ASCII files ``MODEL.mod``)

* ``SliceFile`` class (1d Castro output of the form ``NAME.slice``)

* ``PltFile`` class (h5 files ``pltXXXXX.h5``)

* ``SpectrumFile`` class (h5 files ``SPECTRUM.h5``, produced by sedona if ``output_write_radiation = 1`` is set in the ``param.lua``.)

* ``PlotterClass`` (Visualising DataFiles of all kind either in a set of plots, one combined plot or an animation of the time evolution)

There are also a number of scripts

* ``lightcurve_tools.py`` (An useful tool to create lightcurves from spectrum.h5 files. See section on lightcurve_tools_. for more information eg. how to use it)

* ``test_\*.py`` (test scripts, all available in the ``python_tools/tests/`` folder. They show the basic utilisation of the Datafile/Plotter classes with the example data found in the ``python_tools/data/ directory``).
To use any of the ``test_*.py`` files just type::

  python test_REST_OF_NAME.py

in the ``python_tools/tests/`` folder.

-------------------
Installation
-------------------

To be able to use the tools, make sure, that you have downloaded the python_tools directory to your machine and add the following lines to the start of your python script::

  import sys
  path_to_tools = '../Classes/'
  sys.path.append(path_to_tools)

  import DataFile as DF
  import Plotter as P

``path_to_tools`` can either be a relative path to the python_tools/Classes/ directory or its absolute path (ie. ``/usr/path/to/the/folder/python_tools/Classes/``).

-------------------
File Handling
-------------------
.. _FileHandling:

Generally speaking, Spectrum, Model, Plt and Slice files all have something in common: They consist of a large set of data table and a set of column names - they are some sort of data file.
Examples on how to use them can be found in the python_tools/tests/ directory. There we can see, how to use the classes: ::

  import DataFile as DF

  MyFile = DF.ModelFile(NAME_OF_THE_FILE, PATH_TO_THE_FILE, autoload = False)
  MyFile.load_data()

The ``load_data()`` is not necessary, if ``autoload`` is set to ``True``, which is its default. ``MyFile`` is now a DataFile object (also named DFO). It comes with a list of capabilities,
for example, once we have a file initialised, we can take a quick look at it: ::

  MyFile.plot_data_1D()

This will open up a standard matplotib window. The console and the labels on the plot tell, what is plotted, if the arguments in plot_data_1D are empty, the program will look for the default keys:


.. list-table:: Default ``x`` and ``y`` keys
        :header-rows: 1
        :widths: 20,10,40

        * - File Type
          - ``default_x_key``
          - ``default_y_key``
        * - ModelFile
          - ``'r'``
          - ``'rho'``
        * - SliceFile
          - ``'x'``
          - ``'density'``
        * - PltFile
          - ``'r'``
          - ``'rho'``
        * - SpectrumFile
          - ``'nu'``
          - ``Lnu_averaged'``, which is the automatically angle and time averaged Lnu

The following functions are defined for all data file types:

* ``set_keyword``: Set a value in the data dictionary. ::

    MyFile.set_keyword('numbers',             np.array([1,2,3,4,5]))
    MyFile.set_keyword('other numbers',       np.array([5,4,3,2,1]))
    MyFile.set_keyword('2D array of numbers', np.array([  [1,2,3,4,5],
                                                          [6,7,8,9,10],
                                                          [11,12,13,14,15],
                                                          [16,17,18,19,20],
                                                          [21,22,23,24,25]] ))

.. list-table:: ``set_keyword``
        :header-rows: 1
        :widths: 15,15,40

        * - Argument
          - Type
          - Description
        * - ``key``
          - arbitrary, string recommended
          - key, at which the value is written
        * - ``value``
          - arbitrary, numpy array recommended
          - dictionary entry

* ``get_value``: Return the set value in the data dictionary.::

    numbers = Myfile.get_value('numbers')

.. list-table:: ``get_value``
        :header-rows: 1
        :widths: 15,15,40

        * - Argument
          - Type
          - Description
        * - ``key``
          - arbitrary, string recommended
          - key, at which the value is written

* ``plot_data_1D``: Creates a 1D matplotlib plot of two 1D arrays ``MyFile.data[x]`` vs. ``MyFile.data[y]`` ::

    MyFile.plot_data_1D(x = 'numbers', y = 'other numbers')

* ``plot_data_2D``: Visualises a 2D array vs. two 1D arrays, ie. ``MyFile.data[x]``, ``MyFile.data[y]`` vs. ``MyFile.data[z]`` ::

    MyFile.plot_data_2D(x = 'numbers', y = 'other numbers', z = '2D array of numbers')

.. list-table:: ``plot_data_1D and plot_data_2D``
        :header-rows: 1
        :widths: 15,15,40

        * - Argument
          - Type
          - Description
        * - ``x``, ``y``, ``z`` (for 2D only)
          - arbitrary, string recommended
          - key of the ``x``/``y``/``z`` data in the data dictionary, if no set it will default to the ``default_x/y/z_key``
        * - ``plot_type``
          - string either ``'std'``, ``'logx'``, ``'logy'`` or ``'loglog'``
          - Describes the axis scaling, standard (linear on both axes), semilogarithmic plot with logarithmic x/y axis or loglog plot
        * - ``plot_title``, ``plot_xlabel``, ``plot_ylabel``
          - string
          - Set these by hand to manually write title and labels of the plot. They **do** **not** affect the data plotted!
        * - ``plot_fmt``
          - string
          - format string used by ``matplotib.pyplot.plot()``, eg. ``'co'`` for cyan dots or ``'r--'`` for a dashed red line.
        * - ``plot_save``
          - boolean
          - whether to save the plot or not. Default is ``False``
        * - ``plot_save_type``
          - string
          - What file the plot is to be saved to, eg. ``'pdf'`` or ``'JPG'``. Default is ``'png'``
        * - ``interactive``
          - boolean
          - Whether or not to open up an interface, where two sets of buttons allow to change x and y. See Plotting_ for more information.


But there are also file type specific functions, which are

.. list-table:: File Specific Functions
        :header-rows: 1
        :widths: 20,30,30,30

        * - File Type
          - Function
          - Parameters
          - Description
        * - ``ModelFile``
          - ``add_composition``
          - ``Z``: list of all Z. Default is ``[1]``.
            ``A``: list of all A. Default is ``[1]``.
            ``starting_composition``: list of composition at lowest ``r``. Default is ``[1.]``.
            ``size``: how many entries of the composition have to be calculated.
            ``constant``: Whether to assume constant composition. Default is ``True``.
            ``distrbution_function``: If ``constant`` is set to ``False``, then ``distrbution_function`` (r (float) -> composition (list of floats) ) will be used to calculate the composition at all ``r``.
          - add a composition to the model file
        * -
          - ``get_integrated_mass``
          - n/a
          - integrate rho dr over all available r to find total enclosed mass.
        * -
          - ``write_file_1d``
          - ``type``: Type of output file, either ``'h5'`` or ``'ascii'``. Default is 'ascii',
            ``outpath``: Set, to manually determine, where the file is supposed to be written to.
          - write the information saved in the ``ModelFile`` object to a file usable for sedona.
        * - ``SliceFile``
          - ``find_time``
          - n/a
          - read in the time (float) from the header
        * -
          - ``create_1dmodel_file``
          - ``path_to_model_file``: Use to manually set the path of the ModelFile object
            ``name_of_model_file``: Use to manually set the name of the ModelFile object
            ``size_model_file``: How many data entries the ModelFile will have. The function rebins the original data into the new shape. Default is ``100``.
            ``r_min_model_file``: Manually set the ``r_min`` in the ModelFile. Default is ``0``
            ``v_model_file``: How to determine the velocity profile, either ``"homolog expan"`` for v = r/t or ``"from file"`` for v = p/rho
            ``rebin_model_file``: How to rebin the data. So far only ``'std'`` is available. Other methods are not stable enough
            ``add_composition``: Whether or not to create a pure HI composition to the problem. Default is ``False``.
          - write the information saved in the ``SliceFile`` object into a ModelFile object.
        * - ``PltFile``
          - n/a
          - n/a
          -
        * - ``SpectrumFile``
          - ``find_velocity_coordinates``
          - ``'nu_0'``. Default is Ly alpha frequency.
          - Transforming frequency coord. to velocity coord. via (Delta v)/c (here v/c) = Delta nu/nu_0
        * -
          - ``fit_BB_temperature``
          - ``L``
            ``Tmin = 1e-10``
            ``Tmax = 1e10``
            ``plot_fit = False``
          - Fitting the data (nus, averaged Lnus) to a blackbody spectrum. T is seen as a parameter between Tmin and Tmax. The fitting needs L, which is either given or obtained by using get_integrated_luminosity()
        * -
          - ``subtract_BB_continuum``
          - ``T = None``
            ``L = None``
            ``Tmax = 1e10``
            ``plot_fit = False``
          - Subtract a BB spectrum of temperature T [K] and luminosity L [erg/s] from the continuum
          result is saved as data['Lnu_noBBcont']
        * -
          - ``get_bolometric_lum``
          - n/a
          - returns Lbol (time) ie. Lnu integrated dnu, averaged over viewing angle, still time dependent
        * -
          - ``get_bolometric_mag``
          - n/a
          - returns the bolometric light curve (abs. magnitude)
        * -
          - ``get_ABMag``
          - ``band``
          - returns the AB magnitude light curve, for a given band

        * -
          - ``get_color``
          - ``band1``, ``band2``
          - takes two bands and returns the color
        * -
          - ``get_lcs``
          - ``bands``
          - takes a list of bands and returns their AB magnitude light curves
        * -
          - ``write_bands``
          - ``bands``
          - writes out the file to "lightcurve.out"

-------------------
Tools
-------------------
.. _lightcurve_tools:

In the python_tools/tools directory there is a script called ``lightcurve_tools.py``.
It can be used to create the light curve in a given band.
A list of all available bands can be found in the ``FILTER_LIST`` file in ``python_tools/data/``.
They can also be accessed by typing ::

  python lightcurve_tools.py --bands

To produce a light curve just use: ::

  python lightcurve_tools.py -s <path_to_spectrum_file> -b <space separated list of bands>

The resulting lighcurve can be found in the file ``./lightcurve.out``.
It has 3 columns  called *Time (Days)*, *Lbol (erg/s)*, *Mbol* and one for every band
For example the ``lightcurve.out`` created by running ::

  python lightcurve_tools.py -s ../data/spectrum_1D_TypeIa.h5 -b U B V

in the ``python_tools/tools/`` directory produces a ``lightcurve.out`` (in ``python_tools/tools/``)
with 6 data columns for time, bolometric luminosity, bolometric magnitude and the magnitude in the U, B and V band. All magnitudes are given in the AB Magnitude System.

**Warning:** Invalid elements are assigned a value of 0, the script also floors the magnitudes to 0 if <0

-------------------
Plotting
-------------------
.. _Plotting:

There are also possibilities to compare single or multiple data files. The important basic examples can be found in ``/tests/test_plotter.py`` and ``tests/test_datafilelasses.py``.
The different types of plots available are:

For a single DataFile

* plot X vs. Y, where X and Y are two arrays found in the data

* plot X vs. Y, where two row of buttons let you switch between all available X and Y arrays

* plot X and Y vs. Z

For multiple DataFiles

* making subplots for a set of files (plotting their respective default X vs. default Y)

* combining a set of files into one single plot (again, plotting their respective default X vs. default Y)

* timeseries of a set of files, ie. plotting time dependent data into a short repeating movie. This works for single DataFiles, that have a time axis (ie SpectrumFiles) or a set of files, that have a time associated to them (so far only SliceFiles have been implemented)

A single DataFile plot can be opened without explicitly calling the plotter, for example by::

  path = "python_tools/data/"
  name = "plt19400.slice"

  S = DF.SliceFile(name, path)
  S.plot_data_1D()

The function ``DATAFILE.plot_data_1D()`` takes the arguments described in the FileHandling_ section.
If a 2D plot is needed, use ``DATAFILE.plot_data_2D()``, which can be used to visualise a 2D array of shape (A x B), using two 1D arrays of length A and B as x and y-axis.
For an interactive plotting interface (for 1D plots), set ``interactive = True`` for the ``DATAFILE.plot_data_1D()`` function, ie.: ::

  S.plot_data_1D(interactive = True)

The same plot will be opened, if ::

  Plot = P.Plotter(DataFileObjects = [S])
  Plot.start_1D_UI(mode = "interactiveplot")

is called. The first line initialises the Plotter object called Plot, the second starts the interface.
**Note, that P.Plotter takes a list of DFOs to initialise!** The reason is, that normally one would want to compare a set of DFOs.
Generally speaking, the Plotter class comes (aside from some convenience methods) with a one-for-all-function **``start_1D_UI()``**.
Its most important argument is ``mode``. ``mode`` can be one of the following options:

* ``mode = 'interativeplot'``: Opens the same plot, as if ``DATAFILE.plot_data_1D(interactive)`` were called. Window consists of one plot and two columns of (radio) buttons, that can be used to select the data on the x and y axis.

* ``mode = 'subplots'``: Open a window with N = (number of DFOs in the DataFileObjects list) subplots, one for every DFO. Plots default x vs. default y for every file. There is a button to en- and disable log scaling for the axes.

* ``mode = 'onewindow'``: Open a window with one plot. As for ``mode = 'subplots'`` default x vs. default y is plotted for every file, but now in one combined plot. There is a button to en- and disable log scaling for the axes. **Warning:** The Plotter doesn't check the dimensions and units of the data, it assumes, that it makes sense, that all the data can be put into one plot window!

* ``mode = 'timeseries'``: Open a window with one plot, a progress bar above it, and two buttons for the log scaling. If the Plotter object has been initialised with ``DataFileObjects = [ONEFILE]``, then this file needs to have a time axis. So far, only ``SPECTRUM.h5`` files do. If the ``DataFileObjects`` list consists of more than one element, then it uses ``DFO.find_time()`` to see, if a time can be found.
  Once the data is loaded, the movie is created using the given fps number and film_duration and interpolating between snapshots created from the DFO(s).
