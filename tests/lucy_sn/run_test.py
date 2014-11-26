import os
import pylab as py

def run(pdf,plotup):

    ###########################################
    # run the code
    ###########################################
    os.system("./gomc param_1D.lua")
       
    ###########################################
    # compare the output
    ###########################################
    py.clf()

    x,y,c = py.loadtxt('optical_spectrum.dat',unpack=1,skiprows=1)
    x = x/3600.0/24.0
    py.plot(x,y,'o',color='black',markersize=8,markeredgewidth=2,markerfacecolor='none')

    x,y,c = py.loadtxt('gamma_spectrum.dat',unpack=1,skiprows=1)
    x = x/3600.0/24.0
    py.plot(x,y,'-',color='black')

    # benchmark results
    x,y = py.loadtxt('lucy_lc.dat',unpack=1)
    py.plot(x,y,color='red',linewidth=2)
    x,y = py.loadtxt('lucy_gr.dat',unpack=1)
    py.plot(x,y,color='red',linewidth=2)


    py.title
    py.legend(['sedona LC','sedona GR','lucy LC','lucy GR'])
    pdf.savefig()
    py.xlim(0,60)

    if (plotup):
        py.ion()
        py.show()
        j = raw_input()


    # clean results
    os.system("rm spectrum* gamma* ray_* gomc")
