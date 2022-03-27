# Average surface roughness
# Last modified: June 15 2018
# Coded by Takuro TOKUNAGA

import math
import numpy as np
import cmath
import time
import sys
import random
from scipy.integrate import trapz, simps, quad, quadrature, romberg

start = time.time()

# unit conversion
ucnano = 1.0*np.power(10.,-9)

# plate information
nmax = 20 # average number
nxmax = 100 # use int number
nymax = 100
lxmin = 0
lx = lxmin
lxmax = 700*ucnano
dlx = (lxmax-lxmin)/nxmax
lymin = 0
ly = lymin
lymax = 200*ucnano
dly = (lymax-lymin)/nymax

# parameters
average_heightx = 0

# roughness
height_min = -2*ucnano # [m]
height_max = 2*ucnano # [m]

# file open
f1 = open('average height.txt', 'w') # write mode
f2 = open('roughness1d.txt', 'w') # write mode


# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

for i in range(0, nmax): # average number
    for j in range(0, nxmax+1): # x direction
        # random number
        height = random.uniform(height_min,height_max)

        # Output with text files
        f2.write(str(lx/ucnano)) # x [nm]
        f2.write(str(' '))
        f2.write(str(height/ucnano)) # roughness [nm]
        f2.write('\n')

        # x updatae
        lx = lx  + dlx

        # for average calculation
        average_heightx = average_heightx + height

    # average calculation
    average_heightx = average_heightx/(nxmax+1)

    f1.write(str(i)) # [-]
    f1.write(str(' '))
    f1.write(str(average_heightx/ucnano)) # average_heightx [nm]
    f1.write('\n')

    # reset parameters
    average_heightx = 0
    lx = 0

# file close
f1.close()

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
