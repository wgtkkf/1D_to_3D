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

# number of files
M=25
N=500+M
Np=N-M # N prime
dxmin = -0.5 # [nm]
dx = dxmin # [nm]
dxmax = 1 # [nm]
ddx = (dxmax-dxmin)/(2*N) # [nm]

# parameters
ave = 0.0 # average
stdev = 3.0 # [nm], standard deviation
cl = 17 # correlation length, [nm]

# main
begin()

# file open
f1 = open('correlation_weight.txt', 'w') # write mode
f2 = open('correlation_random.txt', 'w') # write mode

# random number & weight generation
rand = np.zeros(2*N,dtype=np.float64)
weight = np.zeros(2*N,dtype=np.float64)

for i in range(0, 2*N): # M<<N
    rand[i] = normal(ave,stdev) # random number generation, [-] but [nm]
    weight[i] = np.exp(2*np.power(dx/cl,2.0))
    #dx = dx + ddx # [m]

    # Output with text files
    f2.write(str(i)) # [-]
    f2.write(str(' '))
    f2.write(str(rand[i])) # [-], but [nm]
    f2.write(str(' '))
    f2.write(str(weight[i])) # [-], but [nm]
    f2.write('\n')


# re-initialization
shiftx = dxmax*(M/N) # [nm]
dxmin = dxmin + shiftx
dxmax = dxmax - shiftx
ddx = (dxmax-dxmin)/(2*Np)
dx = dxmin

sx = dxmin
dsx = ddx
for i in range(0, 2*Np): # M<<N
    height = 0
    for j in range(0, 2*M):

        height = height + weight[j]*rand[j+i]

    # gaussian correlaiton
    cx = np.exp(-np.power(sx/cl,2.0))

    # Output with text files
    f1.write(str(dx)) # [nm]
    f1.write(str(' '))
    f1.write(str(height)) # [nm]
    f1.write(str(' '))
    f1.write(str(cx)) # [-]
    f1.write('\n')

    dx = dx + ddx # [m]
    sx = sx + dsx

# ends
end()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
