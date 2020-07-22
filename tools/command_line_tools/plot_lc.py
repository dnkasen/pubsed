# Script to test out bokeh plotting sedona spectrum files
#
# To run Type
#    bokeh serve --show  ./plot_spectra.py
# (can add --port 5006 or other number if wanted)

import os
from os.path import dirname, join
import sys

from bokeh.io import curdoc
from bokeh.layouts import row,column
from bokeh.models import ColumnDataSource, Select
from bokeh.plotting import figure
from bokeh.models import Slider
from bokeh.models import Button
from bokeh.models import Range1d
from bokeh.models.widgets import RangeSlider
from bokeh.models.formatters import PrintfTickFormatter
from bokeh.io import output_file, show
from bokeh.models import FileInput

from bokeh.io import output_file, show
from bokeh.layouts import widgetbox
from bokeh.models.widgets import CheckboxButtonGroup

import sedonalib as sed


#




fname = sys.argv[1]
print sys.argv

app_dir = dirname(__file__)

# initial default values
init_time = 20.0
init_lrange = [1000,10000]

bands = ['bol','B','V','R','I','J','H','L','Z']
bands_checkbox = CheckboxButtonGroup(labels=bands, active=[0, 1])

# read spectral data using sedonalib
specdata = sed.SpectrumFile(fname,spec_units='angstrom')

##xx,yy = specdata.get_spectrum(time=init_time)
#full_lrange = [min(xx),max(xx)]
#full_frange = [min(yy),max(yy)]
t,lc= specdata.get_lightcurve()
full_trange = [min(t),max(t)]
full_frange = [min(lc),max(lc)]

print full_trange
# define a time slider for spectra
#time_slider = Slider(start=full_trange[0], end=full_trange[1], value=init_time, step=1, title="time")

# define sliders for setting plot ranges
xrange_slider = RangeSlider(start=full_trange[0], end=full_trange[1], value=(full_trange[0],full_trange[1]), step=100, title="time range")
# define sliders for setting plot ranges
yrange_slider = RangeSlider(start=full_frange[0], end=full_frange[1], value=(full_frange[0],full_frange[1]), step=100, title="y range")
#tick_format = Slider(name='Distance', format=PrintfTickFormatter(format='%.3f m'))

# define the source data container
source = ColumnDataSource(data=dict(x=[], y=[]))

# make the plot
p = figure(plot_height=600, plot_width=800, title="", toolbar_location="below")
p.line(x="x", y="y", source=source, line_width=2)
p.background_fill_color = "#efefef"
#p.x_range=Range1d(1000,9000)
#p.y_range=Range1d(0,50)


# function that updates data selected
def select_data():
    #data_val = data_select.value

#    print plot_type

    #if (data_val != plot_type):
    #    if (data_val == "light curve"):
    #        inputs = column(data_select)
    #    else:
    #        inputs = #column(data_select,time_slider)
    #    curdoc().clear()
    #    curdoc().add_root(row(inputs, p, #width=1100))
    #plot_type = data_val

#    time_val = time_slider.value
    xx,yy = specdata.get_lightcurve()

    return xx,yy

def update_range(newval,type):
    if (type == 1):
        oldval = xrange_slider.value
        xrange_slider.value = (newval,oldval[1])
    if (type == 2):
        oldval = xrange_slider.value
        xrange_slider.value = (oldval[0],newval)

# function to update plot
def update():
    x, y = select_data()
    source.data = dict(x=x, y=y)

    p.x_range.start = xrange_slider.value[0]
    p.x_range.end   = xrange_slider.value[1]
    p.y_range.start = yrange_slider.value[0]
    p.y_range.end   = yrange_slider.value[1]

# define what sliders will do
xrange_slider.on_change('value',lambda attr, old, new: update())
yrange_slider.on_change('value',lambda attr, old, new: update())

p.x_range.on_change('start', lambda attr, old, new: update_range(new,1))
p.x_range.on_change('end',   lambda attr, old, new: update_range(new,2))
p.y_range.on_change('start', lambda attr, old, new: update_range(new,3))
p.y_range.on_change('end',   lambda attr, old, new: update_range(new,4))

inputs = column(xrange_slider,yrange_slider,bands_checkbox)

update()

curdoc().add_root(row(inputs, p, width=1100))
curdoc().title = "Sedona Spectrum"
