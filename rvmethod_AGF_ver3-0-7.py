# Read surface roughness data + liner interporation
# Lego block shape bumpy
# Averaging method
# Last modified: October 08 2018
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
ucpico = 1.0*np.power(10.,-12) # pico to m

# plate information
f0 = open('nxmax.txt', 'r') # read mode (input, fixed value)
for line1 in f0: # read nxmax
    #print(str(line1))
    nxmax = float(line1)

nymax = nxmax
tip_diagonal = 210*ucnano # [m], 170~250 [nm]
tip_length = tip_diagonal/np.sqrt(2) # [m], 120.2~176.7 [nm]
Atip = np.power(tip_length,2) # total area of tip
dAtip = Atip/((nxmax+1)*(nymax+1)) # discretized tip area
total_area = 0 # initialization
area = dAtip

# parameters for atoms
raddi_si = 111*ucpico # silicon [m]
raddi_pt = 177*ucpico # platinum [m]
raddi_ptsi = raddi_si + raddi_pt #  silicon & platinum [m]
lc_si = 5.43*ucangs # [m], Silicon lattice constant
lc_pt = 3.92*ucangs # [m], Platinum lattice constant
lc_ptsi = 0.5*(lc_si+lc_pt) # [m]

# paramters
cutoff = 0.5431 # [nm]
cutoff = lc_ptsi/ucnano # [nm]
bmcounter = 0 # bumpy counter

# gap paramters
gapmin = 0.01 # [nm]
gapmax = 101 # [nm]
gap = gapmin # [nm]
dgap = gapmin # [nm], initial value

# file open
f1 = open('ave.txt', 'r') # read mode
f2 = open('std.txt', 'r') # read mode
f3 = open('conductance_AGF.txt',"w") # write mode
f4 = open('area.txt',"w") # write mode
f5 = open('bumpy_info.txt',"w") # write mode
f6 = open('cond_ratio_AGF.txt',"w") # write mode

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

# read parameters for rand function
for line1 in f1:
    print(str(line1))
    para_ave = float(line1)

for line2 in f2:
    print(str(line2))
    para_std = float(line2)

# read bumpy data
bumpy_2d = pd.read_csv("../surface/bumpy/data0.txt", sep=" ", header=None)
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
    ac_htc = 0
    nc_htc = 0
    total_htc = 0
    ac_conductance = 0
    nc_conductance = 0

    for i in range(0, row): # loop for each bumpy
        bumpy = bumpy_1d[i] # table

        # individual bumpy's gap
        #igap = round(gap-(bumpy+para_std*0.5),5) # [nm]
        igap = round(gap-bumpy,5) # [nm]

        if igap <= cutoff: # [nm]
            igap = cutoff
            bmcounter = bmcounter + 1

            # atomistic contact conductance
            ac_htc = ac_htc + liner_interpolate(igap) # [W/m2K]
        else:
            nc_htc = nc_htc + liner_interpolate(igap) # [W/m2K]

        total_htc = total_htc + liner_interpolate(igap) # [W/m2K]
        total_area = total_area + area

    # output result
    total_conductance = (total_htc*total_area/(((nxmax+1)*(nymax+1))*ucnano)) # [nW/K] ([n*W/m2K]) & averaging
    ac_conductance = (ac_htc*total_area/(((nxmax+1)*(nymax+1))*ucnano)) # [nW/K] ([n*W/m2K]) & averaging
    nc_conductance = (nc_htc*total_area/(((nxmax+1)*(nymax+1))*ucnano)) # [nW/K] ([n*W/m2K]) & averaging

    # ratio of conductance (different from counting bmcounter itself)
    ratio_ac_conductance = (ac_conductance/total_conductance)*100
    ratio_nc_conductance = (nc_conductance/total_conductance)*100

    f3.write(str(gap)) # [nm]
    f3.write(str(' '))
    f3.write(str(total_conductance)) # [nW/K]
    f3.write('\n')

    f4.write(str(gap)) # [nm]
    f4.write(str(' '))
    f4.write(str(total_area)) # [m2]
    f4.write('\n')

    f6.write(str(gap)) # [nm]
    f6.write(str(' '))
    f6.write(str(round(ratio_ac_conductance,5))) # [%]
    f6.write(str(' '))
    f6.write(str(round(ratio_nc_conductance,5))) # [%]
    f6.write(str(' '))
    f6.write(str(100*bmcounter/((nxmax+1)*(nymax+1)))) # [%], contact ratio
    f6.write(str(' '))
    f6.write(str(100*(1-bmcounter/((nxmax+1)*(nymax+1))))) # [%], contact ratio
    f6.write('\n')

    # dgap update
    if gap < 0.098: # ~ 0.1 nm
        dgap = 0.01
    elif gap>0.099 and gap < 0.98: # 0.1 nm ~ 1.0 nm
        dgap = 0.1
    elif gap>0.99 and gap < 9.98: # 1.0 nm ~ 10.0 nm
        dgap = 0.1
    elif gap>=9.99: # 10.0 nm ~
        dgap = 10

    # gap update
    gap = gap + dgap
    gap = round(gap,3)

    # reset parameters
    total_conductance = 0
    ac_htc = 0
    nc_htc = 0
    total_htc = 0
    ac_conductance = 0
    nc_conductance = 0
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
f0.close()
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
