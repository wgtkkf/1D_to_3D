# Surface roughness generator + liner interporation for FE
# Lego block shape bumpy
# Last modified: October 04 2018
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
from lip_FE import liner_interpolate_FE

# Unit conversion:
ucev = 1.602176620898*np.power(10.,-19)
ucnano = 1.0*np.power(10.,-9)
ucangs = 1.0*np.power(10.,-10)

# plate information
f0 = open('nxmax.txt', 'r') # read mode (input, fixed value)
for line1 in f0: # read nxmax
    #print(str(line1))
    nxmax = int(line1)

nymax = nxmax
tip_diagonal = 210*ucnano # [m], 170~250 [nm]
tip_length = tip_diagonal/np.sqrt(2) # [m]
Atip = np.power(tip_length,2) # total area of tip
dAtip = Atip/((nxmax+1)*(nymax+1)) # discretized tip area
total_area = 0 # initialization
area = dAtip

# paramters
cutoff = 0.5431 # [nm]

# gap paramters
gapmin = 0.01 # [nm]
gapmax = 101 # [nm]
gap = gapmin # [nm]
dgap = gapmin # [nm], initial value

temp_low = 300 # [K], low temperature
temp_high = 470 # [K], high temperature
temp_dif = (temp_high)-(temp_low) # [K]

# file open
f1 = open('conductance_FE.txt',"w") # write mode

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

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

    for i in range(0, row): # loop for each bumpy
        bumpy = bumpy_1d[i] # table

        #
        igap = gap-bumpy # [nm], individual bumpy's gap

        if igap < cutoff: # [nm]
            igap = cutoff

        total_conductance = total_conductance + area*liner_interpolate_FE(igap) # [m2]*[W/m2]
        total_area = total_area + area

    # output result
    total_conductance = (total_conductance/ucnano)/temp_dif # [nW/K] from [W]

    f1.write(str(round(gap,2))) # [nm]
    f1.write(str(' '))
    f1.write(str(total_conductance)) # [nW/K]
    f1.write('\n')

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

    # reset parameters
    total_conductance = 0
    total_area = 0
    igap = 0
    bumpy = 0

# file close
f0.close()
f1.close()


finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
