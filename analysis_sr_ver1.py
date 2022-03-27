# Surface roughness analysis
# Last modified: July 08 2019
# Coded by Takuro TOKUNAGA

import numpy as np
import time
import sys
import pandas as pd
import statistics
start = time.time()

nmax = 30
counter = 0 # file number counter
counter2 = 0 #
bumpy_1d_max = np.zeros(nmax,dtype=np.float64)

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

for i in range(0, nmax): # file number
    # heater side
    bumpy_2d = pd.read_csv("../surface/bumpy_tip/rename/data"+str(counter)+".txt", sep=" ", header=None)
    #bumpy_2d = pd.read_csv("../surface/bumpy/data"+str(counter)+".txt", sep=" ", header=None)
    bumpy_2d.columns = ["x", "y", "h"]
    row, col = bumpy_2d.shape # row & column of matorix
    bumpy_1d = np.zeros(row, dtype='float64') # 1024 data points
    bumpy_1d_all = np.zeros(row*nmax, dtype='float64') # all the data points

    for j in range(0, row): # data points loop
        bumpy_1d_all[counter2] = bumpy_2d.iat[j,2] # x line
        counter2 = counter2 + 1

        bumpy_1d[j] = bumpy_2d.iat[j,2] # x line

    # max value
    bumpy_1d_max[i] = np.amax(bumpy_1d)

    # file name counter update
    counter = counter + 1
    #print(str(bumpy_1d_max[i]))

# analysis
stdev = np.std(bumpy_1d_all)
max = np.amax(bumpy_1d_all)
min = np.amin(bumpy_1d_all)

# analysis within max values
max_average = np.average(bumpy_1d_max)
max_std = np.std(bumpy_1d_max)

# display
#print("standard deviatio:{:.2f}".format(stdev) + "[-]")
#print("max:{:.2f}".format(max) + "[nm]")
#print("min:{:.2f}".format(min) + "[nm]")
print("ave of maxes:{:.2f}".format(max_average) + "[nm]")
print("stdev of maxes:{:.2f}".format(max_std) + "[nm]")

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
