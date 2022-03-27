# Average & Std Calculator for rvmethod_series.py
# Coded by Takuro TOKUNAGA
# Last modified: August 29 2018

import numpy as np
import math
import pandas as pd
import matplotlib               # 1. matplotlib for mac
matplotlib.use('TkAgg')         # 2. matplotlib for mac
import matplotlib.pyplot as plt # 3. matplotlib for mac
import matplotlib.font_manager as font_manager
from scipy.stats import norm
from scipy.stats import t # t distribution
from numpy.random import *
import pandas as pd
import time
import statistics
from scipy.stats import gaussian_kde
from scipy.stats import weibull_min
from scipy.stats import exponweib
import seaborn as sns
import pyper
r = pyper.R(use_pandas = "True")
start = time.time()

# function: begin
def begin():
    print ("begin")

# function: end
def end():
    print ("end")

# unit conversion
ucnano = 1.0*np.power(10.,-9)

# main
begin()

# file open
f1 = open('ave.txt', 'r') # read mode
f2 = open('std.txt', 'r') # read mode

# read parameters for rand function
for line1 in f1:
    print(str(line1))
    para_ave = float(line1)

for line2 in f2:
    print(str(line2))
    para_std = float(line2)

## read bumpy data
bumpy_2d = pd.read_csv("../surface/bumpy/data0.txt", sep=" ", header=None)
bumpy_2d.columns = ["x", "y", "h"]
row, col = bumpy_2d.shape # row & column of matorix
bumpy_1d = np.zeros((row), dtype='float64')

# input data into tables
for i in range(0, row):
    bumpy_1d[i] = bumpy_2d.iat[i,2] # x line
    #print(str(bumpy_1d[i]))

# statistical values
min_value_bumpy = min(bumpy_1d) # min [nm]
max_value_zsub = max(bumpy_1d) # max [nm]
mean_value = statistics.mean(bumpy_1d) # mean [nm]
stdev_value = statistics.stdev(bumpy_1d) # stdev [-]

## fitting parameters
param = norm.fit(bumpy_1d) # fitting paramter (average, standarad deviation)
print(param)

# graph
csfont = {'fontname':'Times New Roman'} # define font
plt.figure
x = np.linspace(min_value_bumpy,abs(min_value_bumpy), 100) # range of x axis for fitting curve
pdf_AFM = norm.pdf(x,loc=para_ave, scale=para_std) # gaussian fitting with parameters
pdf_bumpy = norm.pdf(x,loc=param[0], scale=param[1]) # gaussian fitting with parameters
plt.plot(x, pdf_AFM, color = 'red',label="Gaussian (AFM raw data)") # display gaussian fitting with parameters
plt.plot(x, pdf_bumpy, color = 'fuchsia', label="Gaussian (Bumpy data)") # display gaussian fitting with parameters
plt.hist(bumpy_1d, normed=True, color = 'blue', label = 'Bumpy data (Python default)') # histgram of original data, default setting

## graph information
plt.title('Distribution', **csfont) # graph title
plt.xlabel('Bumpy [nm]', fontdict=None, labelpad=None, **csfont)
plt.ylabel('Normalized density [-]', fontdict=None, labelpad=None, **csfont)

# font for legend
font = font_manager.FontProperties(family='Times New Roman',
                                   weight='bold',
                                   style='normal', size=10)
plt.legend(loc='upper right', prop=font) # legend

# plot options
plt.xticks([-5.0,-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0], **csfont)
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], **csfont)

# graph save & display
plt.savefig("histgram2.png") # 1. file saving (1. should be before 2.)
plt.show()                   # 2. file showing (2. should be after 1.)

# file close
f1.close()
f2.close()

end()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
