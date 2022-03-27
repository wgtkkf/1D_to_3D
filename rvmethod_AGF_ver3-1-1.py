# Surface roughness generator + liner interporation
# Rodolphe's method, single results
# Last modified: August 25 2018
# stably running
# Coded by Takuro TOKUNAGA

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
vfactor = (voltage1/voltage2)

tip_length = lxmax
Atip = np.power(tip_length,2) # total area of tip
dAtip = Atip/((nxmax+1)*(nymax+1)) # discretized tip area
division = 51 # division number of each discretized tip area, odd number
dsAtip = dAtip/(division*division) # x & y

cutoff = 0.54 # [nm]
total_area = 0 # initialization
sub_total_area = 0 # initialization

l1 = 0.5*(lxmax/nxmax)
l2 = l1/(division-1)

gapmin = 0.01 # [nm]
gapmax = 101 # [nm]
gap = gapmin # [nm]
dgap = gapmin # [nm], initial value

# 1d tensor data
vx=np.zeros((nxmax*nymax), dtype='float64')
vxcounter = 0

# bumpy counter
bmcounter = 0

# file open
f1 = open('average height.txt', 'w') # write mode
f2 = open('bumpy.txt',"w") # write mode
f3 = open('ave.txt', 'r') # read mode
f4 = open('std.txt', 'r') # read mode
f5 = open('conductance.txt',"w") # write mode
f6 = open('area.txt',"w") # write mode

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
            # random roughness
            bumpy = vfactor*normal(para_ave,para_std) # average & standard deviation

            # Output with text files
            f2.write(str(lx/ucnano)) # x [nm]
            f2.write(str(' '))
            f2.write(str(ly/ucnano)) # y [nm]
            f2.write(str(' '))
            f2.write(str(bumpy)) # roughness [nm]
            f2.write('\n')

            # whole data points
            if vxcounter<(nxmax*nymax):
                vx[vxcounter] = bumpy
            vxcounter = vxcounter+1

            # contact counter
            if (gap-bumpy) < cutoff: # [nm]
                bmcounter = bmcounter + 1

            while sub_total_area < dAtip: # loop for in-plane rotatiaon
                for k in range(0, int((division-1)*0.5)): # each discretized points, 1 direction
                    # height of each discretized area
                    h1 = (l2/l1)*bumpy

                    #
                    area = dsAtip

                    #
                    igap = gap-h1 # [nm]
                    if igap < cutoff:
                        igap = cutoff

                    #
                    total_conductance = total_conductance + area*(liner_interpolate(igap)) # one direction
                    sub_total_area = sub_total_area + area

                    # update of l2
                    l2 = (l1/(division-1))*(k+1)

            #
            total_area = total_area + sub_total_area # close to dAtip

            # y update
            ly = ly  + dly

            # reset parameters
            l2 = l1/(division-1)
            sub_total_area = 0

        # reset parameters
        ly = 0
        l1 = 0.5*(lxmax/nxmax)
        l2 = l1/(division-1)

        # x updatae
        lx = lx  + dlx

    # reset lx & ly
    lx = 0
    ly = 0

    # output result
    total_conductance = (total_conductance/ucnano)

    f5.write(str(gap)) # [nm]
    f5.write(str(' '))
    f5.write(str(total_conductance)) # [nW/K]
    f5.write('\n')

    f6.write(str(gap)) # [nm]
    f6.write(str(' '))
    f6.write(str(total_area)) # [m2]
    f6.write(str(' '))
    f6.write(str(100*bmcounter/((nxmax+1)*(nymax+1)))) # [%]
    f6.write('\n')

    # dgap update
    if gap < 0.098: # ~ 0.1 nm
        dgap = 0.01
    elif gap>0.099 and gap < 0.98: # 0.1 nm ~ 1.0 nm
        dgap = 0.1
    elif gap>0.99 and gap < 10: # 1.0 nm ~ 10.0 nm
        dgap = 1.0
    elif gap>=10: # 10.0 nm ~
        dgap = 10

    # gap update
    gap = gap + dgap

    # reset parameters
    total_conductance = 0
    total_area = 0
    sub_total_area = 0
    igap = 0
    l1 = 0.5*(lxmax/nxmax)
    l2 = l1/(division-1)
    bmcounter = 0

# total average, max
total_mean_value = statistics.mean(vx)
total_max = max(vx)
total_min = min(vx)
f1.write(str(total_mean_value)) # [nm]
f1.write(str(' '))
f1.write(str(total_max)) # [nm]
f1.write(str(' '))
f1.write(str(total_min)) # [nm]
f1.write(str(' '))
f1.write(str(total_max+abs(total_min))) # [nm], peak to peak

# file close
f1.close()
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
