# Surface roughness generator + liner interporation
# Rodolphe's method
# Last modified: August 06 2018
# Coded by Takuro TOKUNAGA
# not working correctly

import math
import numpy as np
import cmath
import time
import sys
import pandas as pd
from scipy.integrate import trapz, simps, quad, quadrature, romberg
from numpy.random import *
import statistics
start = time.time()

# external function
sys.path.append('../regression/')
from lip import liner_interpolate

# Unit conversion:
ucev = 1.602176620898*np.power(10.,-19)
ucnano = 1.0*np.power(10.,-9)
ucangs = 1.0*np.power(10.,-10)

# plate information
nxmax = 32 # use int type
nymax = nxmax

# lx & ly
lxmin = 0
lx = lxmin
lxmax = 200*ucnano
dlx = (lxmax-lxmin)/nxmax
lymin = 0
ly = lymin
lymax = 200*ucnano
dly = (lymax-lymin)/nymax

# parameters
voltage1 = (-28.00) # [nm]
voltage2 = (-30.00) # [nm]
vfactor = voltage1/voltage2

tip_length = lxmax
Atip = np.power(tip_length,2) # total area of tip
dAtip = Atip/((nxmax+1)*(nymax+1)) # discretized tip area
division = 51 # division number of each discretized tip area, odd number
dsAtip = dAtip/division

dangle = (np.pi*0.5)/((division-1)*0.5+1)
angle1 = dangle

gapmin = 0.01 # [nm]
gapmax = 100 # [nm]
gap = gapmin # [nm]
dgap = gapmin # [nm], initial value

temp_low = 300 # [K], low temperature
temp_high = 470 # [K], high temperature
temp_dif = (temp_high)-(temp_low) # [K]

# 1d tensor data
vx=np.zeros((nxmax*nymax), dtype='float64')
vxcounter = 0

# file open
f1 = open('average height.txt', 'w') # write mode
f2 = open('bumpy.txt',"w") # write mode
f3 = open('ave.txt', 'r') # read mode
f4 = open('std.txt', 'r') # read mode
f5 = open('conductance.txt',"w") # write mode

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

# read parameters for rand function
for line3 in f3:
    #print(str(line3))
    para_ave = float(line3)

for line4 in f4:
    #print(str(line4))
    para_std = float(line4)

# distance loop
while gap < gapmax:
    total_conductance = 0
    igap = 0

    for i in range(0, nxmax+1): # x direction
        for j in range(0, nymax+1): # y direction
            # random number
            bumpy = vfactor*normal(para_ave,para_std) # average & standard deviation

            # Output with text files
            f2.write(str(lx/ucnano)) # x [nm]
            f2.write(str(' '))
            f2.write(str(ly/ucnano)) # y [nm]
            f2.write(str(' '))
            f2.write(str(bumpy)) # roughness [nm]
            f2.write('\n')

            if vxcounter<(nxmax*nymax):
                vx[vxcounter] = bumpy
            vxcounter = vxcounter+1

            # individual gap
            igap = bumpy + gap

            # each discretized area
            for i in range(0, int((division-1)*0.5)):
                # height of each discretized area, sphere type
                h1 = igap*np.sin(angle1) # change here

                total_conductance = total_conductance + 2*liner_interpolate(h1)

                # update of l2
                angle1 = angle1 + dangle
                #print(str(angle1))

            # conductance calculation
            total_conductance = total_conductance + liner_interpolate(igap)

            # y update
            ly = ly  + dly

        # reset ly
        ly = 0
        angle1 = dangle

        # x updatae
        lx = lx  + dlx

    # reset lx & ly
    lx = 0
    ly = 0

    # output result
    total_conductance = dsAtip*(total_conductance/ucnano)/temp_dif

    f5.write(str(gap)) # [nm]
    f5.write(str(' '))
    f5.write(str(total_conductance)) # [nW/K]
    f5.write('\n')

    # dgap update
    if gap < 0.1:
        dgap = 0.01
    elif gap>=0.1 and gap < 1.0:
        dgap = 0.1
    elif gap>=1.0 and gap < 10:
        dgap = 1.0
    elif gap>=10:
        dgap = 10

    # gap update
    gap = gap + dgap

    # reset parameters
    total_conductance = 0
    igap = 0
    angle1 = dangle

# total average
total_mean_value = statistics.mean(vx)
f1.write(str(total_mean_value)) # [nm]

# file close
f1.close()
f2.close()
f3.close()
f4.close()
f5.close()

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
