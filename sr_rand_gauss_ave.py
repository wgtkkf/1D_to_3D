# Surface roughness generator
# Last modified: March 12 2018
# Coded by Takuro TOKUNAGA
# Latest version: used for Pt-Si system

import math
import numpy as np
import cmath
import time
import sys
import pandas as pd
from numpy.random import *
import statistics
start = time.time()

# unit conversion
ucnano = 1.0*np.power(10.,-9) # nm to m

# number of files
nmax = 33

# plate information
f0 = open('nxmax.txt', 'r') # read mode (input, fixed value)
for line1 in f0: # read nxmax
    #print(str(line1))
    nxmax = int(line1)

nymax = nxmax

# lx & ly
tip_diagonal = 210*ucnano # [m], 170~250 [nm]
lxmin = 0
lx = lxmin
lxmax = tip_diagonal/np.sqrt(2) # [m], 120.2~176.7 [nm]
dlx = (lxmax-lxmin)/nxmax
lymin = 0
ly = lymin
lymax = lxmax
dly = (lymax-lymin)/nymax

# parameters for surface roughness
counter = 0
overcounter = 0 # confidence interval counter
edged_counter = 0

# 1d tensor data
vx=np.zeros((nxmax*nymax)*nmax, dtype='float64')
vxcounter = 0
voc=np.zeros(nmax, dtype='float64') # vector for overcounter
rms = 0 # rms value

# file open
f1 = open('all_info.txt', 'w') # write mode
# f2: opened below
f3 = open('ave.txt', 'r') # read mode
f4 = open('std.txt', 'r') # read mode

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

random = random()

# read parameters for rand function
for line3 in f3:
    print(str(line3))
    para_ave = float(line3)

for line4 in f4:
    print(str(line4))
    para_std = float(line4)

for i in range(0, nmax):
    f2 = open("../surface/bumpy/data"+str(counter)+".txt","w")
    for j in range(0, nxmax+1): # x direction
        for k in range(0, nymax+1): # y direction
            # gaussian random number
            bumpy = normal(para_ave,para_std) # average & standard deviation

            # sigma manipulation
            cs = 2.17 # sigma coefficient

            if bumpy >=cs*para_std or bumpy <=-cs*para_std:
                overcounter = overcounter + 1

                while True:
                    bumpy = normal(para_ave,para_std)

                    # counter check
                    if bumpy >=cs*para_std or bumpy <=-cs*para_std:
                        overcounter = overcounter + 1

                    if bumpy < cs*para_std or bumpy > -cs*para_std:
                        break

            # bumpy forced substitution
            if bumpy > 4.4:
                bumpy = para_ave
                edged_counter = edged_counter + 1

            # Output with text files
            f2.write(str(lx/ucnano)) # x [nm]
            f2.write(str(' '))
            f2.write(str(ly/ucnano)) # y [nm]
            f2.write(str(' '))
            f2.write(str(bumpy)) # roughness [nm]
            f2.write('\n')

            if vxcounter<(nxmax*nymax)*nmax:
                vx[vxcounter] = bumpy

            ly = ly  + dly
            vxcounter = vxcounter+1

        # reset parameters
        ly = 0
        # x updatae
        lx = lx  + dlx

    # store overcounter to vector
    voc[i] = overcounter

    # display over counter values
    print(str(overcounter))

    # reset parameters
    lx = 0
    ly = 0
    overcounter = 0

    # file name counter update
    counter = counter + 1

# total average
total_mean_height = statistics.mean(vx) # [nm]
total_peak_to_peak = max(vx)-min(vx) # [nm]
total_stdev_height = statistics.stdev(vx) # [-]

# Total RMS value
vx2 = np.square(vx) #
total_rms = np.sum(vx2)/(((nxmax+1)*(nymax+1))*nmax)
total_rms = np.sqrt(total_rms)

# confidence interval
ci = (1-max(voc)/((nxmax+1)*(nymax+1)))*100

# display results
print("Overall mean height:{:.4f}".format(total_mean_height) + "[nm]")
print("Overall peak to peak:{:.3f}".format(total_peak_to_peak) + "[nm]")
print("Overall standard deviation:{:.3f}".format(total_stdev_height) + "[-]")
print("Overall RMS:{:.3f}".format(total_rms) + "[-]")
print("Overall RMS:{:.3f}".format(total_rms) + "[-]")
print("Overall confidence interval:{:.3f}".format(ci) + "[%]")

# file close
f0.close()
f1.close()
f2.close()
f3.close()
f4.close()

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
