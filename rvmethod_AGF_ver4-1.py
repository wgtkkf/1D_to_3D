# Read many surface roughness data files + liner interporation
# Lego block shape bumpy
# Last modified: March 12 2019
# Coded by Takuro TOKUNAGA
# rvmethod, latest code

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
ucnano = 1.0*np.power(10.,-9) # nm to m
ucangs = 1.0*np.power(10.,-10) # angs to m
ucpico = 1.0*np.power(10.,-12) # pico to m

nmax = 33 # file number
counter = 0 # file number counter

# plate information
f0 = open('nxmax.txt', 'r') # read mode (input, fixed value)
for line1 in f0: # read nxmax
    #print(str(line1))
    nxmax = int(line1)

nymax = nxmax

# plate information
tip_diagonal = 210*ucnano # [m], 170~250 [nm]
tip_length = tip_diagonal/np.sqrt(2) # [m]
Atip = np.power(tip_length,2) # total area of tip
dAtip = Atip/((nxmax+1)*(nymax+1)) # discretized tip area
contact_area = 0 # initialization

# parameters for atoms
raddi_si = 111*ucpico # silicon [m]
raddi_pt = 177*ucpico # platinum [m]
raddi_ptsi = 0.5*(raddi_si + raddi_pt) #  silicon & platinum [m]
lc_si = 5.43*ucangs # [m], Silicon lattice constant
lc_pt = 3.92*ucangs # [m], Platinum lattice constant
lc_ptsi = 0.5*(lc_si+lc_pt) # [m]
lc_au = 4.0782*ucangs # [m], Gold lattice constant

# paramters
cutoff = lc_ptsi/ucnano # [nm], pt-si
#cutoff = lc_au/ucnano # [nm], gold
bmcounter = 0 # bumpy counter

# gap paramters
gapmin = 0.01 # [nm]
gapmax = 101 # [nm]
gap = gapmin # [nm]
dgap = gapmin # [nm], initial value

# file open
f1 = open('ave.txt', 'r') # read mode
f2 = open('std.txt', 'r') # read mode

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

# read parameters for rand function
for line1 in f1:
    #print(str(line1))
    para_ave = float(line1)

for line2 in f2:
    #print(str(line2))
    para_std = float(line2)

for l in range(0, nmax): # file number
    bumpy_2d = pd.read_csv("../surface/bumpy/data"+str(counter)+".txt", sep=" ", header=None)
    bumpy_2d.columns = ["x", "y", "h"]
    row, col = bumpy_2d.shape # row & column of matorix
    bumpy_1d = np.zeros(row, dtype='float64')

    # input data into tables
    for j in range(0, row):
        bumpy_1d[j] = bumpy_2d.iat[j,2] # x line
        #print(str(bumpy_1d[i]))

    f3 = open("../surface/conductance/data"+str(counter)+".txt","w")
    f4 = open("../surface/area/data"+str(counter)+".txt","w")
    f5 = open("../surface/ratio/data"+str(counter)+".txt","w")
    f6 = open("../surface/cond_nc_ratio/data"+str(counter)+".txt","w")
    f7 = open("../surface/cond_ac_ratio/data"+str(counter)+".txt","w")

    # distance loop
    while gap < gapmax:
        igap = 0
        total_conductance = 0
        ac_conductance = 0
        nc_conductance = 0

        for i in range(0, row): # loop for each bumpy
            bumpy = bumpy_1d[i] # table

            # individual bumpy's gap
            igap = round(gap-bumpy,5) # [nm]

            if igap <= cutoff: # [nm]
                igap = cutoff
                bmcounter = bmcounter + 1

                # atomistic contact conductance
                ac_conductance = ac_conductance + dAtip*liner_interpolate(igap) # [W/K] ([W/m2K]*[m2])
            else:
                nc_conductance = nc_conductance + dAtip*liner_interpolate(igap) # [W/K] ([W/m2K]*[m2])

            total_conductance = total_conductance + dAtip*liner_interpolate(igap) # [W/K] ([W/m2K]*[m2])
            contact_area = dAtip*bmcounter

        # output result
        total_conductance = (total_conductance/ucnano) # [nW/K] ([n*W/K])
        ac_conductance = (ac_conductance/ucnano) # [nW/K] ([n*W/K])
        nc_conductance = (nc_conductance/ucnano) # [nW/K] ([n*W/K])
        ac_conductance_ratio = ac_conductance/total_conductance*100 # [%]
        nc_conductance_ratio = nc_conductance/total_conductance*100 # [%]
        contact_area_ratio = contact_area/Atip*100 # [%]

        if l==0: # total conductance
            f3.write(str(gap)) # [nm]
            f3.write(str(' '))
            f3.write(str(total_conductance)) # [nW/K]
            f3.write('\n')
        else:
            f3.write(str(total_conductance)) # [nW/K]
            f3.write('\n')

        if l==0: # contacted area ratio
            f4.write(str(gap)) # [nm]
            f4.write(str(' '))
            f4.write(str(contact_area_ratio)) # [%]
            f4.write('\n')
        else:
            f4.write(str(contact_area_ratio)) # [%]
            f4.write('\n')

        if l==0: # contact ratio
            f5.write(str(gap)) # [nm]
            f5.write(str(' '))
            f5.write(str(100*bmcounter/((nxmax+1)*(nymax+1)))) # [%]
            f5.write('\n')
        else:
            f5.write(str(100*bmcounter/((nxmax+1)*(nymax+1)))) # [%]
            f5.write('\n')

        if l==0: # non-contact ratio of conductance
            f6.write(str(gap)) # [nm]
            f6.write(str(' '))
            f6.write(str(nc_conductance_ratio)) # [%]
            f6.write('\n')
        else:
            f6.write(str(nc_conductance_ratio)) # [%]
            f6.write('\n')

        if l==0: # atomistic-contact ratio of conductance
            f7.write(str(gap)) # [nm]
            f7.write(str(' '))
            f7.write(str(ac_conductance_ratio)) # [%]
            f7.write('\n')
        else:
            f7.write(str(ac_conductance_ratio)) # [%]
            f7.write('\n')

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
        nc_conductance = 0
        contact_area = 0
        igap = 0
        bumpy = 0
        bmcounter = 0

    # reset all the parameters for next file
    gap = gapmin
    total_conductance = 0
    ac_conductance = 0
    nc_conductance = 0
    contact_area = 0
    igap = 0
    bumpy = 0
    bmcounter = 0

    # file name counter update
    counter = counter + 1

# file close
f0.close()
f1.close()
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()
f7.close()

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
