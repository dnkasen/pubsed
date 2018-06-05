import pylab as py
import sys
import h5py 


def compare(pdf):


    ###########################################
    # compare the output
    ###########################################
    py.clf()

    fin = h5py.File('optical_spectrum.h5')
    t = py.array(fin['time'])
    L = py.array(fin['Lnu'])
    mu = py.array(fin['mu'])
    t = t/3600.0/24.0
    for i in range(len(mu)):
      py.plot(t,L[:,0,i])

    # gamma-ray deposition
 #   tdep = []
 #   gdep = []
 #   for j in range(1,150):
 #       ray = 'ray_00'
 #       if (j< 100): ray = ray + '0' 
 #       if (j < 10): ray = ray + '0' 
 #       ray = ray + str(j) + '.dat'
#        fin = open(ray,'r')
 #       line = fin.readline()
  #      tdep.append(float(line.split()[3])/3600.0/24.0)
   #     sum = 0
    #    r0 = 0 
     #   line = fin.readline()
      #  for line in fin:
       #     data = line.split()
        #    r1 = float(data[0])
         #   vol = 4.0*3.14159/3.0*(r1**3 - r0**3)
          #  sum += float(data[5])*vol
        #    r0 = r1
        #gdep.append(sum)
        
#    py.plot(tdep,gdep,'o')
        
    # benchmark results
    x,y = py.loadtxt('lucy_lc.dat',unpack=1)
    py.plot(x,y,color='red',linewidth=2)
    x,y = py.loadtxt('lucy_gr.dat',unpack=1)
    py.plot(x,y,color='blue',linewidth=2)
    py.ylim(1e40,0.4e44)

    py.title
    py.legend(['sedona LC','sedona GR','lucy LC','lucy GR'])
    py.xlim(0,60)

    if (pdf != ''): pdf.savefig()
    else:
        py.ion()
        py.show()
        j = raw_input()

if __name__=='__main__': compare('')

