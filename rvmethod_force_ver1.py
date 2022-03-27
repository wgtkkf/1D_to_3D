# Read many surface roughness data files + external function for force model
# Lego block shape bumpy
# Last modified: March 08 2019
# Coded by Takuro TOKUNAGA
# rvmethod, latest code

import math
import numpy as np
import cmath
import time
import sys
import pandas as pd
start = time.time()

# external function for force model
sys.path.append('../force/')
from lj_frc import dlj
from lifshitz_frc_asym import lifshitz_force_asym
from cap_plate_frc import cap_force
from es_frc import coulomb_force

# Unit conversion:
ucev = 1.602176620898*np.power(10.,-19) # electron volt to joule
ucnano = 1.0*np.power(10.,-9) # nano to m
ucangs = 1.0*np.power(10.,-10) # angs to m
ucpico = 1.0*np.power(10.,-12) # pico to m

# parameters for surface roughness file
nmax = 33 # file number
counter = 0 # file number counter

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
contact_area = 0 # initialization

# parameters for atoms
raddi_si = 111*ucpico # silicon [m]
raddi_pt = 177*ucpico # platinum [m]
raddi_ptsi = 0.5*(raddi_si + raddi_pt) #  silicon & platinum [m]
Aptsi = np.pi*np.power(raddi_ptsi,2.0) # [m2]
lc_si = 5.43*ucangs # [m], Silicon lattice constant
lc_pt = 3.92*ucangs # [m], Platinum lattice constant
lc_ptsi = 0.5*(lc_si+lc_pt) # [m]
lc_au = 4.0782*ucangs # [m], Gold lattice constant

# paramters
cutoff = lc_ptsi/ucnano # [nm]
#cutoff = lc_au/ucnano # [nm]

# gap paramters
gapmin = 5.0 # [nm]
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

# read parameters for rand function of surface roughness (not in this code)
for line1 in f1:
    #print(str(line1))
    para_ave = float(line1)

for line2 in f2:
    #print(str(line2))
    para_std = float(line2)

for l in range(0, nmax): # file number
    # reading generated surface roughness
    bumpy_2d = pd.read_csv("../surface/bumpy/data"+str(counter)+".txt", sep=" ", header=None)
    bumpy_2d.columns = ["x", "y", "h"]
    row, col = bumpy_2d.shape # row & column of matorix
    bumpy_1d = np.zeros(row, dtype='float64')

    # input data into tables
    for j in range(0, row):
        bumpy_1d[j] = bumpy_2d.iat[j,2] # x line

    f3 = open("../surface/force_total/data"+str(counter)+".txt","w") # total force
    f4 = open("../surface/force_lj/data"+str(counter)+".txt","w") # Lennard-Jones force
    f5 = open("../surface/force_casimir/data"+str(counter)+".txt","w") # Casimir force
    f6 = open("../surface/force_capacitor/data"+str(counter)+".txt","w") # Capacitor force
    f7 = open("../surface/force_coulomb/data"+str(counter)+".txt","w") # Coulomb force

    # gap distance loop
    while gap < gapmax:
        igap = 0
        lj_f = 0
        casimir_f = 0
        capacitor_f = 0
        coulomb_f = 0
        total_f = 0

        for i in range(0, row): # loop for each bumpy
            bumpy = bumpy_1d[i] # table

            # individual bumpy's gap
            igap = round(gap-bumpy,5) # [nm]
            igap = igap*ucnano # [m], unit conversion

            # call force model
            lj_f =  lj_f + dlj(igap)*dAtip/Atip                     # [N] ([N]*[1/m2]*[m2])
            casimir_f = casimir_f + lifshitz_force_asym(igap)*dAtip # [N] ([N/m2]*[m2])
            capacitor_f = capacitor_f + cap_force(igap)             # [N]
            coulomb_f = coulomb_f + coulomb_force(igap)             # [N]

        total_f = lj_f + casimir_f + capacitor_f + coulomb_f # [N]

        # output result
        lj_f = (lj_f/ucnano) # [nN]
        casimir_f = (casimir_f/ucnano) # [nN]
        capacitor_f = (capacitor_f/ucnano) # [nN]
        coulomb_f = (coulomb_f/ucnano) # [nN]
        total_f = (total_f/ucnano) # [nN]

        if l==0: # total force
            f3.write(str(gap)) # [nm]
            f3.write(str(' '))
            f3.write(str(total_f)) # [nN]
            f3.write('\n')
        else:
            f3.write(str(total_f)) # [nN]
            f3.write('\n')

        if l==0: # Lennard-Jones force
            f4.write(str(gap)) # [nm]
            f4.write(str(' '))
            f4.write(str(lj_f)) # [nN]
            f4.write('\n')
        else:
            f4.write(str(lj_f)) # [nN]
            f4.write('\n')

        if l==0: # Casimir force
            f5.write(str(gap)) # [nm]
            f5.write(str(' '))
            f5.write(str(casimir_f)) # [nN]
            f5.write('\n')
        else:
            f5.write(str(casimir_f)) # [nN]
            f5.write('\n')

        if l==0: # Capacitor force
            f6.write(str(gap)) # [nm]
            f6.write(str(' '))
            f6.write(str(capacitor_f)) # [nN]
            f6.write('\n')
        else:
            f6.write(str(capacitor_f)) # [nN]
            f6.write('\n')

        if l==0: # Coulomb force
            f7.write(str(gap)) # [nm]
            f7.write(str(' '))
            f7.write(str(coulomb_f)) # [nN]
            f7.write('\n')
        else:
            f7.write(str(coulomb_f)) # [nN]
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

        # reset quantities
        lj_f = 0 # [nN]
        casimir_f = 0 # [nN]
        capacitor_f = 0 # [nN]
        coulomb_f = 0 # [nN]
        total_f = 0 # [nN]

    # reset all the parameters for next file
    gap = gapmin # do not forget this line
    igap = 0
    bumpy = 0
    
    lj_f = 0 # [nN]
    casimir_f = 0 # [nN]
    capacitor_f = 0 # [nN]
    coulomb_f = 0 # [nN]
    total_f = 0 # [nN]

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
