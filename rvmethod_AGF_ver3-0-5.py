# Read surface roughness data + liner interporation
# Lego block shape bumpy
# Last modified: September 04 2018
# Coded by Takuro TOKUNAGA

import math
import numpy as np
import cmath
import time
import sys
import pandas as pd
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

# factors
factor_area = 1.5

# plate information
nxmax = 32 # use int type
nymax = nxmax
tip_length = 200*ucnano # [m]
Atip = np.power(tip_length,2) # total area of tip
dAtip = Atip/((nxmax+1)*(nymax+1)) # discretized tip area
total_area = 0 # initialization

# paramters
cutoff = 0.54 # [nm]
bmcounter = 0 # bumpy counter

# gap paramters
gapmin = 0.01 # [nm]
gapmax = 101 # [nm]
gap = gapmin # [nm]
dgap = gapmin # [nm], initial value

# file open
f1 = open('ave.txt', 'r') # read mode
f2 = open('std.txt', 'r') # read mode
f3 = open('conductance.txt',"w") # write mode
f4 = open('area.txt',"w") # write mode
f5 = open('bumpy_info.txt',"w") # write mode

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

# read parameters for rand function
for line1 in f1:
    #print(str(line3))
    para_ave = float(line1)

for line2 in f2:
    #print(str(line4))
    para_std = float(line2)

# read bumpy data
bumpy_2d = pd.read_csv("../surface/bumpy/data1.txt", sep=" ", header=None)
bumpy_2d.columns = ["x", "y", "h"]
row, col = bumpy_2d.shape # row & column of matorix
bumpy_1d = np.zeros((row), dtype='float64')

# input data into tables
for i in range(0, row):
    bumpy_1d[i] = bumpy_2d.iat[i,2] # x line
    #print(str(bumpy_1d[i]))

# distance loop
while gap < gapmax:
    igap = 0
    total_conductance = 0

    for i in range(0, row): # loop for each bumpy
        bumpy = bumpy_1d[i] # table

        # contact counter
        if (gap-bumpy) < cutoff: # [nm]
            bmcounter = bmcounter + 1

        #
        area = dAtip

        #
        igap = gap-bumpy # [nm], individual bumpy's gap

        if igap < cutoff: # [nm]
            igap = cutoff
            area = dAtip*factor_area

        total_conductance = total_conductance + area*liner_interpolate(igap) # [W/K] ([W/m2K]*[m2])
        total_area = total_area + area

    # output result
    total_conductance = (total_conductance/ucnano) # [nW/K] ([n*W/K])

    f3.write(str(gap)) # [nm]
    f3.write(str(' '))
    f3.write(str(total_conductance)) # [nW/K]
    f3.write('\n')

    f4.write(str(gap)) # [nm]
    f4.write(str(' '))
    f4.write(str(total_area)) # [m2]
    f4.write(str(' '))
    f4.write(str(100*bmcounter/((nxmax+1)*(nymax+1)))) # [%]
    f4.write('\n')

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
    igap = 0
    bumpy = 0
    bmcounter = 0

# bumpy information
bumpy_average = statistics.mean(bumpy_1d) # average
bumpy_max = max(bumpy_1d) # max
bumpy_min = min(bumpy_1d) # min
bumpy_pp = bumpy_max+abs(bumpy_min) # peak to peak
bumpy_stdev = statistics.stdev(bumpy_1d) # stdev
# rms
bumpy_rms = 0
for i in range(0, row):
    bumpy_rms = bumpy_rms + np.power(bumpy_1d[i], 2)
bumpy_rms = np.sqrt(bumpy_rms/row)

f5.write(str(bumpy_average)) # [nm]
f5.write('\n')
f5.write(str(bumpy_max)) # [nm]
f5.write('\n')
f5.write(str(bumpy_min)) # [nm]
f5.write('\n')
f5.write(str(bumpy_pp)) # [nm]
f5.write('\n')
f5.write(str(bumpy_stdev)) # [nm]
f5.write('\n')
f5.write(str(bumpy_rms)) # [nm]

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
