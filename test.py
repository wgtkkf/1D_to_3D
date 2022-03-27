# Read many surface roughness data files + liner interporation
# Last modified: September 29 2019
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

def finish():
    print ("finish")

# Unit conversion:
ucnano = 1.0*np.power(10.,-9) # nm to m
ucangs = 1.0*np.power(10.,-10) # angs to m
ucpico = 1.0*np.power(10.,-12) # pico to m

begin()


FE = pd.read_csv("../surface/FE_datasets.txt", sep=" ", header=None)
row, col = FE.shape # row & column of matorix
FE_1d = np.zeros(row, dtype='float64')
for j in range(0, row):
    FE_1d[j] = FE.iat[j,2] # x line
    print(str(FE_1d[j]))

#print(str(row))
#print(str(FE))

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
