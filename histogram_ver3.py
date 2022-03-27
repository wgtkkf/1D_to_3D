# Average & Std Calculator for surface roughness generator (sr_rand_gauss_ave.py)
# use raw data of only afm_snom_8107.txt
# Surface roughness for Pt heater
# Last modified: October 04 2018
# Coded by Takuro TOKUNAGA

import numpy as np
import math
import pandas as pd
import matplotlib               # 1. matplotlib for mac
matplotlib.use('TkAgg')         # 2. matplotlib for mac
import matplotlib.pyplot as plt # 3. matplotlib for mac
import matplotlib.font_manager as font_manager
from scipy.stats import norm
from scipy.stats import skewnorm
from scipy.stats import t # t distribution
from numpy.random import *
import pandas as pd
import time
import statistics
from scipy.stats import gaussian_kde
from scipy.stats import weibull_min
from scipy.stats import exponweib
import scipy.stats as stats
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

voltage1_8107 = (-28.00) # [nm]
voltage2_8107 = (-30.00) # [nm]
vfactor_8107 = (voltage1_8107/voltage2_8107)

# main
begin()

# file open
f1 = open('ave.txt', 'w') # write mode (outputs)
f2 = open('std.txt', 'w') # write mode (outputs)

## read raw data table, AFM_8107
data_2d = pd.read_csv("afm_snom_8107_21.txt", sep="\t", header=0, index_col=0)
#data_2d = pd.read_csv("../surface/AFM_data/AFM_SNOM_8368_SR_Data_python.txt", sep="\t", header=0, index_col=0)
row, col = data_2d.shape # row & column of matorix
data_1d_8107 = np.ravel(data_2d) # 2d to 1d

# statistical values
mean_value = statistics.mean(data_1d_8107) # mean

# paramters
zshift = mean_value
data_1d_zsub_8107=np.zeros((row*col), dtype='float64') # sub: z shift subtract

# data manipulation: 0 point shift & unit conversion [nm]
for i in range(0,row*col):
    data_1d_zsub_8107[i] = vfactor_8107*(data_1d_8107[i] - zshift)/ucnano # convert AFM raw data to Height data

## merge 7825 & 8107
data_1d_zsub = data_1d_zsub_8107

max_value_zsub = max(data_1d_zsub) # max, maybe [nm]
min_value_zsub = min(data_1d_zsub) # min, maybe [nm]
mean_value_zsub = statistics.mean(data_1d_zsub) # mean, maybe [nm]
stdev_value_zsub = statistics.stdev(data_1d_zsub) # stdev

#print(str(max_value_zsub))
#print(str(min_value_zsub))
#print(str(mean_value_zsub))
#print(str(stdev_value_zsub))

## fitting parameters
param = norm.fit(data_1d_zsub) # fitting paramter (average, standarad deviation)
# paramters output
f1.write(str(param[0])) # average
f2.write(str(param[1])) # standard deviation

### for numerical bumpy start ###
# file open
f3 = open('ave.txt', 'r') # read mode
f4 = open('std.txt', 'r') # read mode

### for origin graph ###
f5 = open('origin_histogram.txt', 'w') # read mode
f6 = open('origin_norm.txt', 'w') # read mode

# read parameters for rand function
for line1 in f3:
    print(str(line1))
    para_ave = float(line1)

for line2 in f4:
    print(str(line2))
    para_std = float(line2)

# read numerical bumpy data
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
param_bumpy = norm.fit(bumpy_1d) # fitting paramter (average, standarad deviation)

### for generated bumpy end ###

# graph display
csfont = {'fontname':'Times New Roman'} # define font
plt.figure
x = np.linspace(-6.0,6.0, 100) # range of x axis for fitting curve
pdf_afm = norm.pdf(x,loc=param[0], scale=param[1]) # gaussian fitting with parameters
pdf_bumpy = norm.pdf(x,loc=param_bumpy[0], scale=param_bumpy[1]) # gaussian fitting with parameters
pdf_skewbumpy = skewnorm.pdf(x, a = 0.4, loc=param_bumpy[0], scale=param_bumpy[1]) # gaussian fitting with parameters
plt.plot(x, pdf_afm, 'black', label="Gaussian distribution") # display gaussian fitting with parameters
#plt.plot(x, pdf_bumpy, 'b-', linestyle='dashed', label="Gaussian (Numerical bumpy)") # display gaussian fitting with parameters
#plt.plot(x, pdf_skewbumpy, 'g-', linestyle='dashed', label="Skew Gaussian (Numerical bumpy)") # display gaussian fitting with parameters
plt.hist(data_1d_zsub, normed=True, color = 'red', label = 'AFM surface roughness') # histgram of AFM data, default setting
#plt.hist(bumpy_1d, normed=True, alpha=.7,color = 'blue', label = 'Numerical bumpy') # histgram of original data, default setting

## graph information
#plt.title('Distribution', **csfont) # graph title
plt.xlabel('Height, H [nm]', fontdict=None, labelpad=None, **csfont)
plt.ylabel('Normalized density [-]', fontdict=None, labelpad=None, **csfont)

# font for legend
font = font_manager.FontProperties(family='Times New Roman',
                                   weight='bold',
                                   style='normal', size=10)
plt.legend(loc='upper right', prop=font) # legend

# plot options
plt.xticks([-6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0], **csfont)
plt.yticks([0.0, 0.1, 0.2, 0.3, 0.4], **csfont)
# graph save & display
plt.savefig("histogram.png") # 1. file saving (1. should be before 2.)
plt.show()                  # 2. file showing (2. should be after 1.)

## Q-Q plot
#f = plt.figure(figsize=(12, 8))
#ax = f.add_subplot(111)
#ax.set_title("Probplot")
#stats.probplot(data_1d_zsub, dist='norm', plot=ax) # AFM
#stats.probplot(bumpy_1d, dist='norm', plot=ax) # numerical
#plt.savefig("ppline.png") # 1. file saving (1. should be before 2.)
#plt.show()

#print(stats.shapiro(data_1d_zsub))
#print(stats.shapiro(bumpy_1d))

## output for origin
## output of [0] frequency, [1] class value
ret = plt.hist(data_1d_zsub, normed=True, color = 'red', label = 'AFM surface roughness') # histgram of AFM data, default setting
f5.write(str(ret[0])) # frequency
f5.write('\n')
f5.write('\n')
f5.write(str(ret[1])) # class value

## output for origin
## output of norm
f6.write(str(x)) # x
f6.write('\n')
f6.write('\n')
f6.write(str(pdf_afm)) # norm

# file close
f1.close()
f2.close()
f3.close()
f4.close()
f5.close()

end()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
