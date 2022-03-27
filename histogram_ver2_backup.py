# Average & Std Calculator for rvmethod_series.py
# Coded by Takuro TOKUNAGA
# Last modified: September 28 2018

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
f1 = open('ave.txt', 'w') # write mode (outputs)
f2 = open('std.txt', 'w') # write mode (outputs)

## read data table, AFM_7825
data_2d = pd.read_csv("afm_snom_7825.txt", sep="\t", header=0, index_col=0)
row, col = data_2d.shape # row & column of matorix
data_1d_7825 = np.ravel(data_2d) # 2d to 1d

# statistical values
mean_value = statistics.mean(data_1d_7825) # mean

# paramters
zshift = mean_value
data_1d_zsub_7825=np.zeros((row*col), dtype='float64') # sub: z shift subtract

# data manipulation: 0 point shift & unit conversion [nm]
for i in range(0,row*col):
    data_1d_zsub_7825[i] = (data_1d_7825[i] - zshift)/ucnano

## read data table, AFM_8107
data_2d = pd.read_csv("afm_snom_8107.txt", sep="\t", header=0, index_col=0)
row, col = data_2d.shape # row & column of matorix
data_1d_8107 = np.ravel(data_2d) # 2d to 1d

# statistical values
mean_value = statistics.mean(data_1d_8107) # mean

# paramters
zshift = mean_value
data_1d_zsub_8107=np.zeros((row*col), dtype='float64') # sub: z shift subtract

# data manipulation: 0 point shift & unit conversion [nm]
for i in range(0,row*col):
    data_1d_zsub_8107[i] = (data_1d_8107[i] - zshift)/ucnano

## merge 7825 & 8107
data_1d_zsub = np.append(data_1d_zsub_7825,data_1d_zsub_8107)

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

### for generated bumpy ###
# file open
f3 = open('ave.txt', 'r') # read mode
f4 = open('std.txt', 'r') # read mode

# read parameters for rand function
for line1 in f3:
    para_ave = float(line1)
    
for line2 in f4:
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
param_bumpy = norm.fit(bumpy_1d) # fitting paramter (average, standarad deviation)

### for generated bumpy ###

# graph
csfont = {'fontname':'Times New Roman'} # define font
plt.figure
x = np.linspace(min_value_zsub,abs(min_value_zsub), 100) # range of x axis for fitting curve
#pdf = norm.pdf(x) # gaussian fitting
#plt.plot(x,pdf, 'b-', label="gaussian") # gaussian fitting
pdf_afm = norm.pdf(x,loc=param[0], scale=param[1]) # gaussian fitting with parameters
pdf_bumpy = norm.pdf(x,loc=param_bumpy[0], scale=param_bumpy[1]) # gaussian fitting with parameters
plt.plot(x, pdf_afm, 'r-', label="Gaussian (AFM bumpy)") # display gaussian fitting with parameters
plt.plot(x, pdf_bumpy, 'b-', label="Gaussian (Numerical bumpy)") # display gaussian fitting with parameters
#plt.plot(x, pdf_fitted2, 'g-', label="Gaussian (Self-defined)") # display gaussian fitting with parameters
plt.hist(data_1d_zsub, normed=True, color = 'red', label = 'AFM bumpy (Python default)') # histgram of original data, default setting
plt.hist(bumpy_1d, normed=True, alpha=.7,color = 'blue', label = 'Numerical bumpy (Python default)') # histgram of original data, default setting

## self-defined distribution
#sns.kdeplot(data_1d_zsub, color = 'green', bw=0.75,label="bw: 0.75")
#pdf_fitted2 = norm.pdf(x,loc=param[0], scale=6) # gaussian fitting with parameters
#pdf_weib = exponweib.pdf(x, a=3.5, c=1.2, loc=param[0],scale=param[1]) # gaussian fitting with parameters
#plt.plot(x, pdf_weib, 'g-', label="Weibull (self-defined)") # display gaussian fitting with parameters

# random number generation
rand_weib = exponweib.rvs(a=0.45, c=3.6,loc=0,scale=3.5,size=1000)
#plt.hist(rand_weib-2.0, normed=True, histtype='stepfilled',alpha=.3, color='green', label='Weibull random') # not use density but use normed

## graph information
plt.title('Distribution', **csfont) # graph title
plt.xlabel('AFM data [nm]', fontdict=None, labelpad=None, **csfont)
plt.ylabel('Normalized density [-]', fontdict=None, labelpad=None, **csfont)


# font for legend
font = font_manager.FontProperties(family='Times New Roman',
                                   weight='bold',
                                   style='normal', size=10)
plt.legend(loc='upper right', prop=font) # legend

# plot options
#plt.xticks([-3.5, -2.5, -1.5, 0.0, 1.5, 2.5, 3.5], **csfont)
#plt.xticks([-5.0,-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0], **csfont)
plt.xticks([-6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0], **csfont)
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], **csfont)

# graph save & display
plt.savefig("histgram.png") # 1. file saving (1. should be before 2.)
plt.show()                  # 2. file showing (2. should be after 1.)

# file close
f1.close()
f2.close()
f3.close()
f4.close()

end()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
