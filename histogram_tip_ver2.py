# Average & Std Calculator for surface roughness generator (sr_rand_gauss_ave_tip.py)
# Surface roughness for Si tip
# Last modified: June 17 2019
# Coded by Takuro TOKUNAGA

import numpy as np
import pandas as pd
import matplotlib               # 1. matplotlib for mac
matplotlib.use('TkAgg')         # 2. matplotlib for mac
import matplotlib.pyplot as plt # 3. matplotlib for mac
import matplotlib.font_manager as font_manager
from scipy.stats import norm
from numpy.random import *
import time
import statistics
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

voltage1 = (-28.00) # [nm/V]
voltage2 = (-19.57) # [nm/V]
vfactor = (voltage1/voltage2)

# main
begin()

# file open
f1 = open('ave_tip.txt', 'w') # write mode (outputs)
f2 = open('std_tip.txt', 'w') # write mode (outputs)

## read raw data table, afm_snom_####
data_2d = pd.read_csv("AFM_SNOM_8845_FWD.txt", sep="\t", header=0, index_col=0)
row, col = data_2d.shape # row & column of matorix
data_1d_tip = np.ravel(data_2d) # 2d to 1d

# statistical values
mean_value = statistics.mean(data_1d_tip) # mean

# paramters
zshift = mean_value
data_1d_zsub_temp=np.zeros((row*col), dtype='float64') # sub: z shift subtract

# data manipulation: 0 point shift & unit conversion [nm]
for i in range(0,row*col):
    data_1d_zsub_temp[i] = vfactor*(data_1d_tip[i] - zshift)/ucnano # convert AFM raw data to Height data

#
data_1d_zsub = data_1d_zsub_temp

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
print(str(param[0]))
print(str(param[1]))

### for drawing origin graph ###
f5 = open('origin_histogram_tip.txt', 'w') # write mode
f6 = open('origin_norm_tip.txt', 'w') # write mode

### for numerical bumpy: start ###
# read numerical bumpy data
bumpy_2d = pd.read_csv("../surface/bumpy_tip/data0.txt", sep=" ", header=None)
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

# gaussian fitting parameters for numerical bumpy
param_bumpy = norm.fit(bumpy_1d) # fitting paramter (average, standarad deviation)
#print(str(param_bumpy[0]))
#print(str(param_bumpy[1]))
### for numerical bumpy: end ###

# graph display
csfont = {'fontname':'Times New Roman'} # define font
plt.figure

# for gaussian distribution
x = np.linspace(-1.5,1.5,100) # range of x axis for fitting curve
#pdf_afm = norm.pdf(x,loc=param[0], scale=param[1]) # gaussian fitting with parameters
sf = np.sqrt(2*np.pi*np.power(param[1],2.0)) # scale factor
#print(str(param[1]))
pdf_afm = 0.19*norm.pdf(x,loc=param[0], scale=param[1]) # gaussian fitting with parameters

## for raw data histogram
weights = np.ones_like(data_1d_zsub)/len(data_1d_zsub)
bin_row = np.sqrt(row)

## for numerical data histogram
weights_numerical = np.ones_like(bumpy_1d)/len(bumpy_1d)

plt.plot(x, pdf_afm, color = 'black', label="Gaussian fitting") # display gaussian fitting with parameters
plt.hist(data_1d_zsub, weights=weights,bins=int(bin_row), color = 'red', label = 'Measured roughness') # histgram of AFM data, default setting
plt.hist(bumpy_1d, weights=weights_numerical,bins=10, alpha=0.3, color = 'blue', label = 'Numerical roughness') # histgram of AFM data, default setting

## graph information
#plt.title('Distribution', **csfont) # graph title
plt.xlabel('Tip Surface Roughness, H [nm]', fontdict=None, labelpad=None, **csfont)
plt.ylabel('Normalized density [-]', fontdict=None, labelpad=None, **csfont)

# font for legend
font = font_manager.FontProperties(family='Times New Roman',
                                   weight='bold',
                                   style='normal', size=10)
plt.legend(loc='upper right', prop=font) # legend

# plot options
plt.xticks([-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5], **csfont)
plt.yticks([0.0, 0.1, 0.2, 0.3, 0.4], **csfont)
# graph save & display
plt.savefig("histogram.png") # 1. file saving (1. should be before 2.)
plt.show()                  # 2. file showing (2. should be after 1.)

## output for origin
## output of [0] frequency, [1] class value
ret = plt.hist(data_1d_zsub, normed=True, color = 'red', label = 'Measured roughness') # histgram of AFM data, default setting
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
f5.close()
f6.close()

end()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
