import pylab as py
import sys



def compare(pdf):


    ###########################################
    # compare the output
    ###########################################
    py.clf()

    x,y,c = py.loadtxt('optical_spectrum.dat',unpack=1,skiprows=1)
    x = x/3600.0/24.0
    py.plot(x,y,'o',color='black',markersize=8,markeredgewidth=2,markerfacecolor='none')

#    x,y,c = py.loadtxt('gamma_spectrum.dat',unpack=1,skiprows=1)
#    x = x/3600.0/24.0
#    py.plot(x,y,'-',color='black')

    # gamma-ray deposition
    tdep = []
    gdep = []
    for j in range(1,150):
        ray = 'ray_00'
        if (j< 100): ray = ray + '0'
        if (j < 10): ray = ray + '0'
        ray = ray + str(j)
        fin = open(ray,'r')
        line = fin.readline()
        tdep.append(float(line.split()[3])/3600.0/24.0)
        sum = 0
        r0 = 0 
        line = fin.readline()
        for line in fin:
            data = line.split()
            r1 = float(data[0])
            vol = 4.0*3.14159/3.0*(r1**3 - r0**3)
            sum += float(data[5])*vol
            r0 = r1
            print vol
        gdep.append(sum)
        
    py.plot(tdep,gdep,'o')
        
    # benchmark results
    x,y = py.loadtxt('lucy_lc.dat',unpack=1)
    py.plot(x,y,color='red',linewidth=2)
    x,y = py.loadtxt('lucy_gr.dat',unpack=1)
    py.plot(x,y,color='blue',linewidth=2)


    py.title
    py.legend(['sedona LC','sedona GR','lucy LC','lucy GR'])
    py.xlim(0,60)

    if (pdf != ''): pdf.savefig()
    else:
        py.ion()
        py.show()
        j = raw_input()

compare('')
