# Surface roughness generator for tip side, Kim et al., Nature, 528 387 2015
# Last modified: August 11 2019
# Coded by Takuro TOKUNAGA

import numpy as np
import time
import sys
import pandas as pd
from numpy.random import *
import statistics
start = time.time()

# functions for begin & finish
def begin():
    print ("begin")

def end():
    print ("finish")

# unit conversion
ucnano = 1.0*np.power(10.,-9) # nm to m

# parameters
N=10000 # number of data points
ave = 0.0 # average
#stdev = 15.0 # [nm], standard deviation
stdev = 3.0 # [nm], standard deviation
cl = 17 # correlation length, [nm]

# random number & weight generation
rand = np.zeros(N,dtype=np.float64)

# main
begin()

# file open
f1 = open('../surface/bumpy_gold/data0.txt', 'w') # write mode

for i in range(0, N): # M<<N
    rand[i] = normal(ave,stdev) # random number generation, [-] but [nm]

    # Output with text files
    f1.write(str(i)) # [-], dummy
    f1.write(str(' '))
    f1.write(str(i)) # [-], dummy
    f1.write(str(' '))
    f1.write(str(rand[i])) # [-], but [nm]
    f1.write('\n')

# correlation function origion
czero = np.exp

total_mean_height = statistics.mean(rand) # [nm]
total_max = max(rand)
total_min = min(rand)
total_peak_to_peak = total_max-total_min # [nm]

rand2 = np.square(rand) # square
total_rms = np.sum(rand2)/N # sum & devision
total_rms = np.sqrt(total_rms) # sqrt

print("Overall mean height:{:.4f}".format(total_mean_height) + "[nm]")
print("Overall max height:{:.4f}".format(total_max) + "[nm]")
print("Overall min height:{:.4f}".format(total_min) + "[nm]")
print("Overall peak to peak:{:.3f}".format(total_peak_to_peak) + "[nm]")
print("Overall RMS:{:.3f}".format(total_rms) + "[-]")

# ends
end()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
