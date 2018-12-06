from numpy import pi, sin
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import h5py
import bisect




def get_spec(time,Lnu,tget):
	indt = bisect.bisect(time,float(tget))
	print tget,indt
	return Lnu[indt,:]

# read spectrum data
fname = 'spectrum.h5'
fin = h5py.File(fname,'r')
time = np.array(fin['time'],dtype='d')
Lnu  = np.array(fin['Lnu'],dtype='d')
nu   = np.array(fin['nu'],dtype='d')
time = time/3600.0/24.0
fin.close()

# convert to Llam
L = Lnu*nu**2/2.99e10/1e8
x  = 2.99e10/nu*1e8

t0 = 20.0
freq_0 = 1

xlimits = [1000,10000]

axis_color = 'lightgoldenrodyellow'

fig = plt.figure()
ax = fig.add_subplot(111)

# Adjust the subplots region to leave some space for the sliders and buttons
fig.subplots_adjust(left=0.25, bottom=0.25)


# Draw the initial plot
# The 'line' variable is used for modifying the line later
[line] = ax.plot(x,get_spec(time,L,t0), linewidth=2, color='red')
ax.set_xlim(xlimits)
#ax.set_ylim([-10, 10])

# Add two sliders for tweaking the parameters

# Define an axes area and draw a slider in it
time_slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03], axisbg=axis_color)
time_slider = Slider(time_slider_ax, 'time (days)', min(time), max(time), valinit=t0)

# Draw another slider
freq_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axis_color)
freq_slider = Slider(freq_slider_ax, 'Freq', 0.1, 30.0, valinit=freq_0)


# Define an action for modifying the line when any slider's value changes
def sliders_on_changed(val):
    line.set_ydata(get_spec(time,L,time_slider.val))
    fig.canvas.draw_idle()
time_slider.on_changed(sliders_on_changed)
freq_slider.on_changed(sliders_on_changed)

# Add a button for resetting the parameters
reset_button_ax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
reset_button = Button(reset_button_ax, 'Reset', color=axis_color, hovercolor='0.975')
def reset_button_on_clicked(mouse_event):
    freq_slider.reset()
    time_slider.reset()
reset_button.on_clicked(reset_button_on_clicked)

# Add a set of radio buttons for changing color
color_radios_ax = fig.add_axes([0.025, 0.5, 0.15, 0.15], axisbg=axis_color)
color_radios = RadioButtons(color_radios_ax, ('red', 'blue', 'green'), active=0)
def color_radios_on_clicked(label):
    line.set_color(label)
    fig.canvas.draw_idle()
color_radios.on_clicked(color_radios_on_clicked)

plt.show()
