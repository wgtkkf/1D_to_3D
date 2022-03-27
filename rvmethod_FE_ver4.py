# Surface roughness generator + liner interporation for FE
# Consider both of tip & heater surface roughness
# Lego block shape bumpy
# Last modified: September 30 2019
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
ucpico = 1.0*np.power(10.,-12) # pico to m

nmax = 45 # file number
counter = 0 # file number counter

# plate information
f0 = open('nxmax.txt', 'r') # read mode (input, fixed value)
for line1 in f0: # read nxmax
    #print(str(line1))
    nxmax = int(line1)

nymax = nxmax
tip_diagonal = 240*ucnano # [m], 170~250 [nm]
tip_length = tip_diagonal/np.sqrt(2) # [m]
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
lc_au = 4.0782*ucangs # [m], Gold lattice constant

# paramters
cutoff = lc_ptsi/ucnano # [nm]
#cutoff = lc_au/ucnano # [nm]

# gap paramters
gapmin = 0.01 # [nm]
gapmax = 101 # [nm]
gap = gapmin # [nm]
dgap = gapmin # [nm], initial value

temp_low = 300 # [K], low temperature
temp_high = 470 # [K], high temperature
temp_dif = (temp_high)-(temp_low) # [K]

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

# input data into tables
for l in range(0, nmax): # file number
    # heater side
    bumpy_2d = pd.read_csv("../surface/bumpy/data"+str(counter)+".txt", sep=" ", header=None)
    bumpy_2d.columns = ["x", "y", "h"]
    row, col = bumpy_2d.shape # row & column of matorix
    bumpy_1d = np.zeros(row, dtype='float64')

    # tip side
    bumpy_2d_tip = pd.read_csv("../surface/bumpy_tip/data"+str(counter)+".txt", sep=" ", header=None)
    bumpy_2d_tip.columns = ["x", "y", "h"]
    row_tip, col_tip = bumpy_2d_tip.shape # row & column of matorix
    bumpy_1d_tip = np.zeros(row_tip, dtype='float64')

    if row != row_tip:
        print('Stopped!')
        break

    # input data into tables
    for j in range(0, row):
        bumpy_1d[j] = bumpy_2d.iat[j,2] # x line
        #print(str(bumpy_1d[i]))

    f1 = open("../surface/conductance_FE/data"+str(counter)+".txt","w")

    # distance loop
    while gap < gapmax:
        igap = 0
        total_conductance = 0

        for i in range(0, row): # loop for each bumpy
            bumpy = bumpy_1d[i] # table
            bumpy_tip = bumpy_1d_tip[i] # table, tip side

            # individual bumpy's gap
            igap = round(gap-bumpy-bumpy_tip,5) # [nm], sign change is considered

            if igap < cutoff: # [nm]
                igap = cutoff

            total_conductance = total_conductance + area*liner_interpolate_FE(igap) # [m2]*[W/m2]
            total_area = total_area + area

        # output result
        total_conductance = (total_conductance/ucnano)/temp_dif # [nW/K] from [W]

        if l==0: # total conductance
            f1.write(str(gap)) # [nm]
            f1.write(str(' '))
            f1.write(str(total_conductance)) # [nW/K]
            f1.write('\n')
        else:
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
        gap = round(gap,3)

        # reset parameters
        total_conductance = 0
        total_area = 0
        igap = 0
        bumpy = 0
        bumpy_tip = 0

    # reset all the parameters for next file
    gap = gapmin
    total_conductance = 0
    total_area = 0
    igap = 0
    bumpy = 0
    bumpy_tip = 0

    # file name counter update
    counter = counter + 1

# file close
f0.close()
f1.close()

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
