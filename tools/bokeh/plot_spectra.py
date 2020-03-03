# Script to test out bokeh plotting sedona spectrum files
#
# To run Type
#    bokeh serve --show  ./plot_spectra.py
# (can add --port 5006 or other number if wanted)
import os
from os.path import dirname, join


from bokeh.io import curdoc
from bokeh.layouts import row,column
from bokeh.models import ColumnDataSource, Select
from bokeh.plotting import figure
from bokeh.models import Slider
from bokeh.models import Button

import sedonalib as sed
specdata = sed.SpectrumFile("toyIa_spectrum.h5")


app_dir = dirname(__file__)

time_slider = Slider(start=0.1, end=100, value=20, step=1, title="time")

options = ['spectrum', 'light curve']
plot_type = options[0]

data_select = Select(title="Plot Type:", value=options[0],
                     options=options)

redraw_button = Button(label='ReDraw')

def change_click():
    print('I was clicked')

redraw_button.on_click(change_click)

source = ColumnDataSource(data=dict(x=[], y=[]))

p = figure(plot_height=600, plot_width=800, title="", toolbar_location="below")
p.line(x="x", y="y", source=source, line_width=2)
p.background_fill_color = "#efefef"


def select_data():
    data_val = data_select.value

#    print plot_type

    #if (data_val != plot_type):
    #    if (data_val == "light curve"):
    #        inputs = column(data_select)
    #    else:
    #        inputs = #column(data_select,time_slider)
    #    curdoc().clear()
    #    curdoc().add_root(row(inputs, p, #width=1100))
    #plot_type = data_val

    time_val = time_slider.value
    xx,yy = specdata.get_spectrum(time=time_val)
    return xx,yy


def update():
    x, y = select_data()
    source.data = dict(x=x, y=y)

data_select.on_change('value', lambda attr, old, new: update())
time_slider.on_change('value', lambda attr, old, new: update())



inputs = column(data_select,time_slider,redraw_button)

update()

curdoc().add_root(row(inputs, p, width=1100))
curdoc().title = "Sedona Spectrum"
