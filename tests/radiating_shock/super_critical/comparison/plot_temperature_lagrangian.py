import numpy as np
import matplotlib
from matplotlib.pyplot import *

Tl1_file = open("ray_00000")
Tl2_file = open("ray_00004")
Tl3_file = open("ray_00008")
Tl4_file = open("ray_00012")
Tl5_file = open("ray_00016")
Tl6_file = open("ray_00020")
Tl7_file = open("ray_00024")
Tl8_file = open("ray_00028")

Tl1 = np.loadtxt(Tl1_file)
Tl2 = np.loadtxt(Tl2_file)
Tl3 = np.loadtxt(Tl3_file)
Tl4 = np.loadtxt(Tl4_file)
Tl5 = np.loadtxt(Tl5_file)
Tl6 = np.loadtxt(Tl6_file)
Tl7 = np.loadtxt(Tl7_file)
Tl8 = np.loadtxt(Tl8_file)


rl1 = Tl1[::,0] - 2.4e12 * np.ones(len(Tl1[::,0]))
Tlg1 = Tl1[::,3] 

rl2 = Tl2[::,0] - 2.4e12 * np.ones(len(Tl1[::,0]))
Tlg2 = Tl2[::,3] 

rl3 = Tl3[::,0] - 2.4e12 * np.ones(len(Tl1[::,0]))
Tlg3 = Tl3[::,3] 

rl4 = Tl4[::,0] - 2.4e12 * np.ones(len(Tl1[::,0]))
Tlg4 = Tl4[::,3]

rl5 = Tl5[::,0] - 2.4e12 * np.ones(len(Tl1[::,0]))
Tlg5 = Tl5[::,3]

rl6 = Tl6[::,0] - 2.4e12 * np.ones(len(Tl1[::,0]))
Tlg6 = Tl6[::,3]

rl7= Tl7[::,0] - 2.4e12 * np.ones(len(Tl1[::,0]))
Tlg7 = Tl7[::,3]

rl8= Tl8[::,0] - 2.4e12 * np.ones(len(Tl1[::,0]))
Tlg8 = Tl8[::,3]

plot(rl2,Tlg2, '-',c = 'b')
plot(rl4,Tlg4, '-',c = 'b')
plot(rl6,Tlg6, '-',c = 'b')
plot(rl8,Tlg8, '-',c = 'b')


xlabel(r"radius (cm)",fontsize=16)
xlim(0,6.e10)

ylabel(r"gas temperature (Kelvin)",fontsize=16)
ylim(0,6.e3)

show()
close()


