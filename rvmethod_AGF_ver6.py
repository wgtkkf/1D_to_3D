# Read many surface roughness data files + liner interporation
# Consider both of tip & heater surface roughness
# Lego block shape bumpy
# Last modified: September 29 2019
# Coded by Takuro TOKUNAGA
# steably working

import math
import numpy as np
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
ucnano = 1.0*np.power(10.,-9) # nm to m
ucangs = 1.0*np.power(10.,-10) # angs to m
ucpico = 1.0*np.power(10.,-12) # pico to m

nmax = 45 # file number
counter = 0 # file number counter

# plate information
f0 = open('nxmax.txt', 'r') # read mode (input, fixed value)
for line1 in f0: # read nxmax
    #print(str(line1))
    nxmax = int(line1)

nymax = nxmax

# plate information
tip_diagonal = 240*ucnano # [m], 170~250 [nm]
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
#print(cutoff*10) # [angs]
#cutoff = lc_au/ucnano # [nm], gold
bmcounter = 0 # bumpy counter
fecounter = 0

# gap paramters
gapmin = 0.01 # [nm]
gapmax = 101 # [nm]
gap = gapmin # [nm]
dgap = gapmin # [nm], initial value

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

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

    # FE part
    FE = pd.read_csv("../surface/FE_datasets.txt", sep=" ", header=None)
    row_fe, col_fe = FE.shape # row & column of matorix
    FE_1d = np.zeros(row_fe, dtype='float64')
    for j in range(0, row_fe):
        FE_1d[j] = FE.iat[j,l+1] # [nW/K]

    if row != row_tip:
        print('Stopped!')
        break

    # input data into tables
    for j in range(0, row):
        bumpy_1d[j] = bumpy_2d.iat[j,2] # x line
        bumpy_1d_tip[j] = bumpy_2d_tip.iat[j,2] # x line
        #print(str(bumpy_1d[i]))
        #print(str(bumpy_1d_tip[i]))

    # outputs
    f3 = open("../surface/conductance/data"+str(counter)+".txt","w")
    f4 = open("../surface/area/data"+str(counter)+".txt","w")
    f5 = open("../surface/ratio/data"+str(counter)+".txt","w")
    f6 = open("../surface/cond_nc_ratio/data"+str(counter)+".txt","w")
    f7 = open("../surface/cond_ac_ratio/data"+str(counter)+".txt","w")

    # distance loop
    while gap < gapmax:
        igap = 0
        agf_conductance = 0
        total_conductance = 0
        ac_conductance = 0
        nc_conductance = 0

        for i in range(0, row): # loop for each bumpy
            bumpy = bumpy_1d[i] # table, heater side
            bumpy_tip = bumpy_1d_tip[i] # table, tip side

            # individual bumpy's gap
            igap = round(gap-bumpy-bumpy_tip,5) # [nm], sign change is considered

            if igap <= cutoff: # [nm]
                igap = cutoff
                bmcounter = bmcounter + 1

                # atomistic contact conductance
                ac_conductance = ac_conductance + dAtip*liner_interpolate(igap) # [W/K] ([W/m2K]*[m2])
            else:
                nc_conductance = nc_conductance + dAtip*liner_interpolate(igap) # [W/K] ([W/m2K]*[m2])

            agf_conductance = agf_conductance + dAtip*liner_interpolate(igap) # [W/K] ([W/m2K]*[m2])
            contact_area = dAtip*bmcounter

        # output result
        total_conductance = (agf_conductance/ucnano) + FE_1d[fecounter] # [nW/K] ([n*W/K]), AGF + FE
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
        agf_conductance = 0
        nc_conductance = 0
        contact_area = 0
        igap = 0
        bumpy = 0
        bumpy_tip = 0
        bmcounter = 0
        fecounter = fecounter + 1

    # reset all the parameters for next file
    gap = gapmin
    total_conductance = 0
    agf_conductance = 0
    ac_conductance = 0
    nc_conductance = 0
    contact_area = 0
    igap = 0
    bumpy = 0
    bumpy_tip = 0
    bmcounter = 0
    fecounter = 0

    # file name counter update
    counter = counter + 1

# file close
f0.close() # input
f3.close() # output
f4.close() # output
f5.close() # output
f6.close() # output
f7.close() # output

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
