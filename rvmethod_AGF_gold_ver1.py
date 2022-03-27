# Read many surface roughness data files + liner interporation
# Consider tip/Not consider heater surface roughness
# Lego block shape bumpy & Spherical surface > Currently no consideration of surface roughness
# Last modified: August 15 2019
# Coded by Takuro TOKUNAGA

import numpy as np
import time
import sys
import pandas as pd
import statistics
start = time.time()

# external function
sys.path.append('../regression/')
from lip import liner_interpolate

# Unit conversion:
ucnano = 1.0*np.power(10.,-9) # nm to m
ucangs = 1.0*np.power(10.,-10) # angs to m
ucpico = 1.0*np.power(10.,-12) # pico to m

nmax = 1 # file number
counter = 0 # file number counter

# plate information
f0 = open('nxmax.txt', 'r') # read mode (input, fixed value)
for line1 in f0: # read nxmax
    #print(str(line1))
    nxmax = int(line1)

nxmax = 100
nymax = nxmax

# plate information
tip_radious = 30*ucnano # [m]
Atip = np.power(tip_radious,2) # total area of tip
dAtip = Atip/((nxmax+1)*(nymax+1)) # discretized tip area

# parameters for radius
nrmax = 10000 # division nrmax toward r direction, integer type
srmin = 0 # [m] initialization of rmin, s:small
sr = srmin # [m] initialization of r, s:small
srmax = tip_radious # [m] initialization of rmax, s:small
dsr = (srmax-srmin)/nrmax # [m] initialization of dr, s:small

# parameters for atoms
lc_au = 4.0782*ucangs # [m], Gold lattice constant

# paramters
cutoff = lc_au/ucnano # [nm], gold
bmcounter = 0 # bumpy counter

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

def dtilde(r, d):
    y = d + srmax - np.sqrt(np.power(srmax,2)-np.power(r,2))

    return y

# main
begin()

for l in range(0, nmax): # file number

    # tip side
    bumpy_2d = pd.read_csv("../surface/bumpy_gold/data"+str(counter)+".txt", sep=" ", header=None)
    bumpy_2d.columns = ["x", "y", "h"]
    row, col = bumpy_2d.shape # row & column of matorix
    bumpy_1d = np.zeros(row, dtype='float64')

    # input data into tables
    for j in range(0, row):
        bumpy_1d[j] = bumpy_2d.iat[j,2] # x line
        #print(str(bumpy_1d[i]))

    # peak value
    bumpy_1d_max = max(bumpy_1d)

    f1 = open("../surface/conductance_AGF_gold/data"+str(counter)+".txt","w")

    # distance loop
    while gap < gapmax:
        igap = 0
        total_conductance = 0

        # initialization of the derjaguin parameters
        sr = srmin
        dtldr = 0
        ind_h = 0
        total_conductance = 0
        bumpy_counter = 0

        # loop for r
        while sr < srmax: # [m]
            # d tilder (dtldr)
            dtldr = dtilde(sr, gap*ucnano) # [m]

            # surface roughness
            bumpy = bumpy_1d[bumpy_counter] # [nm], table
            # individual bumpy's gap
            #if bumpy == bumpy_1d_max:
            #    igap = round(gap,5) # [nm], sign change is considered
            #else:
            #    igap = round(gap+(bumpy_1d_max-bumpy),5) # [nm], sign change is considered

            #igap = round(dtldr/ucnano-bumpy,5) # [nm], sign change is considered
            igap = round(dtldr/ucnano,5) # [nm], sign change is considered

            if igap < cutoff: # [nm]
                igap = cutoff # [nm]

            # individual heat transfer coefficient
            ind_h = liner_interpolate(igap) # ([nm]), [W/m2K]

            # individual conductance
            ind_g = ind_h*(2*np.pi*sr)*dsr # [W/K], (W/m2K)*m2

            # conductance (summation)
            total_conductance = total_conductance + ind_g # [W/K]

            # sr update
            sr = sr + dsr # [m]

            bumpy_counter = bumpy_counter + 1 # [-]

        # output result
        total_conductance = (total_conductance/ucnano) # [nW/K] ([n*W/K])

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
        igap = 0
        bumpy = 0
        bumpy_counter = 0

    # reset all the parameters for next file
    gap = gapmin
    total_conductance = 0
    ac_conductance = 0
    nc_conductance = 0
    igap = 0
    bumpy = 0
    bumpy_counter = 0

    # file name counter update
    counter = counter + 1

# file close
f0.close()
f1.close()

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
