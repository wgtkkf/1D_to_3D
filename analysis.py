# Output files analysis code
# Last modified: August 21 2018
# https://qiita.com/segavvy/items/51a515c19bcd29b13b7f
# Coded by Takuro TOKUNAGA

import math
import numpy as np
import cmath
import time
import sys
import pandas as pd
from scipy.integrate import trapz, simps, quad, quadrature, romberg
from numpy.random import *
import statistics
import glob
import json,csv
start = time.time()

# plate information
nmax = 20 # file number
nxmax = 40
nymax = nxmax
nxymax = nxmax*nymax

#
number = 37 # data points, 0 to 36
sx1=np.zeros((number*nxymax), dtype='float64')
sy1=np.zeros((number*nxymax), dtype='float64')
fx1=np.zeros((number*nxymax), dtype='float64')
sx2=np.zeros((number), dtype='float64')
fx2=np.zeros((number), dtype='float64')
sx3=np.zeros((number), dtype='float64')
fx3=np.zeros((number), dtype='float64')


# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

# file manipulation
for i in range(0, nmax): # file loop
    data1 = pd.read_csv("../surface/bumpy/data"+str(i)+".txt", sep=" ", header=None) # bumpy
    data2 = pd.read_csv("../surface/area/data"+str(i)+".txt", sep=" ", header=None) # area
    data3 = pd.read_csv("../surface/conductance/data"+str(i)+".txt", sep=" ", header=None) # conductance

    # data format
    data1.columns = ["x[nm]", "y[nm]", "roughness[nm]"] # bumpy
    data2.columns = ["gap[nm]", "total area[m2]"] # area
    data3.columns = ["gap[nm]", "total conductance[nW/K]"] # conductance

    # file open
    f3 = open("../surface/bumpy/bumpy/data"+str(i)+".csv", 'w') # write mode
    for j in range(0, nxymax):
        sx1[j] = data1.iat[j,0] # x line
        sy1[j] = data1.iat[j,1] # y line
        fx1[j] = data1.iat[j,2] # fx line

        f3.write(str(' '))
        f3.write(str(fx1[j])) # roughness [nm]
        f3.write('\n')

    # file open
    f4 = open("../surface/area/area/data"+str(i)+".csv", 'w') # write mode
    f5 = open("../surface/conductance/conductance/data"+str(i)+".csv", 'w') # write mode

    # input data into tables
    for j in range(0, number):
        sx2[j] = data2.iat[j,0] # x line
        fx2[j] = data2.iat[j,1] # fx line

        sx3[j] = data3.iat[j,0] # x line
        fx3[j] = data3.iat[j,1] # fx line

        f4.write(str(' '))
        f4.write(str(fx2[j])) # area [m2]
        f4.write('\n')

        f5.write(str(' '))
        f5.write(str(fx3[j])) # conductance [nW/K]
        f5.write('\n')

# file merge, bumpy
csv_files = glob.glob("../surface/bumpy/*.csv")
list = []
for f in csv_files:
    list.append(pd.read_csv(f))
df = pd.concat(list, axis=1)
df.to_csv("../surface/bumpy/total.csv")

# file merge, conductance
csv_files = glob.glob("../surface/conductance/*.csv")
list = []
for f in csv_files:
    list.append(pd.read_csv(f))
df = pd.concat(list, axis=1)
df.to_csv("../surface/conductance/total.csv")

# calculation

# file close
f3.close()
f4.close()
f5.close()

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
