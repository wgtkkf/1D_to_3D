# Surface roughness generator for heater side
# Last modified: June 03 2019
# Coded by Takuro TOKUNAGA
# Latest version: used for Pt-Si system

import numpy as np
import time
import sys
import pandas as pd
from numpy.random import *
import statistics
start = time.time()

# unit conversion
ucnano = 1.0*np.power(10.,-9) # nm to m

# number of files
nmax = 5

# plate information
f0 = open('nxmax.txt', 'r') # read mode (input, fixed value)
for line1 in f0: # read nxmax
    #print(str(line1))
    nxmax = int(line1)

nymax = nxmax

# lx & ly
tip_diagonal = 240*ucnano # [m], 170~250 [nm]
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
edge_counter = 0

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
    para_ave = np.float64(line3)
    print("Average:{:.3f}".format(para_ave) + "[nm]")

for line4 in f4:
    para_std = np.float64(line4)
    print("Standard deviation:{:.3f}".format(para_std) + "[nm]")

for i in range(0, nmax):
    f2 = open("../surface/bumpy/data"+str(counter)+".txt","w")
    for j in range(0, nxmax+1): # x direction
        for k in range(0, nymax+1): # y direction

            # sigma manipulation
            cs = 2.80 # sigma coefficient

            while True:
                # gaussian random number
                bumpy = normal(para_ave,para_std)

                # C.I. filter
                if bumpy > -cs*para_std and bumpy < cs*para_std:
                    break
                elif bumpy >=cs*para_std or bumpy <=-cs*para_std:
                    overcounter = overcounter + 1

            # Edge filter
            if bumpy > 4.47: # [nm] #
                bumpy = para_ave
                overcounter = overcounter + 1
                edge_counter = edge_counter + 1

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
total_points = (nxmax+1)*(nymax+1)*nmax
total_rms = np.sum(vx2)/total_points
total_rms = np.sqrt(total_rms)

# average confidence interval
#ci = (1-np.average(voc)/((nxmax+1)*(nymax+1)))*100
# minimum confidence interval
ci = (1-max(voc)/((nxmax+1)*(nymax+1)))*100


# display results
print("Overall mean height:{:.4f}".format(total_mean_height) + "[nm]")
print("Overall peak to peak:{:.3f}".format(total_peak_to_peak) + "[nm]")
print("Overall standard deviation:{:.3f}".format(total_stdev_height) + "[-]")
print("Overall RMS:{:.3f}".format(total_rms) + "[-]")
print("Overall confidence interval:{:.3f}".format(ci) + "[%]")
print("Overall edge counter:{:.3f}".format((edge_counter/total_points)*100) + "[%]")

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
