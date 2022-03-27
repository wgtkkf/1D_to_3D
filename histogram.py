# Average & Std Calculator for sr_rand_gauss_ave.py
# Coded by Takuro TOKUNAGA
# Last modified: July 10 2018

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

f = open('ave.txt', 'w') # write mode (outputs)
f1 = open('std.txt', 'w') # write mode (outputs)

# read data table
data_2d = pd.read_csv("afm_snom_7825.txt", sep="\t", header=0, index_col=0)
row, col = data_2d.shape # row & column of matorix
data_1d = np.ravel(data_2d) # 2d to 1d

## statistical values
max_value = max(data_1d) # max
min_value = min(data_1d) # min
mean_value = statistics.mean(data_1d) # mean
median_valu = np.median(data_1d)
stdev_value = statistics.stdev(data_1d) # standard deviation

# paramters
zshift = mean_value
data_1d_zsub=np.zeros((row*col), dtype='float64') # sub: z shift subtract

# data manipulation: 0 point shift & unit conversion [nm]
for i in range(0,row*col):
    data_1d_zsub[i] = (data_1d[i] - zshift)/ucnano

max_value_zsub = max(data_1d_zsub)
min_value_zsub = min(data_1d_zsub)
mean_value_zsub = statistics.mean(data_1d_zsub) # mean
stdev_value_zsub = statistics.stdev(data_1d_zsub) # stdev

#print(max_value_zsub)
#print(min_value_zsub)
#print(mean_value_zsub)

#bin_number = 10
#bin_width = 0.5

bin_number = int(math.log2(row*col)+1)-1
#print(str(bin_number))
bin_width = (max_value_zsub-min_value_zsub)/bin_number
#print(str(bin_width))

## fitting parameters
param = norm.fit(data_1d_zsub) # fitting paramter (average, standarad deviation)
#print(param)

# graph
csfont = {'fontname':'Times New Roman'} # define font
plt.figure
x = np.linspace(min_value_zsub,abs(min_value_zsub), 100) # range of x axis for fitting curve
pdf = norm.pdf(x) # gaussian fitting
pdf_fitted = norm.pdf(x,loc=param[0], scale=param[1]) # gaussian fitting with parameters

#plt.plot(x,pdf, 'b-', label="gaussian") # display fitting curve
plt.plot(x, pdf_fitted, 'r-', label="Gaussian (parametric)") # display fitting curve with parameters
#plt.hist(data_1d_zsub, normed=True, label = 'AFM data') # histgram of original data
plt.hist(data_1d_zsub, bins=bin_number, rwidth=bin_width, normed=True, color='blue',label = 'AFM data (Sturges formula based)') # histgram of original data
#h = plt.hist(np.random.triangular(-2.5, 0, 2.5, 1000),bins=200,normed=True,color='green',label="Triangle") # triangle

plt.title('Distribution', **csfont) # graph title
plt.xlabel('AFM data [nm]', fontdict=None, labelpad=None, **csfont)
plt.ylabel('Normalized density [-]', fontdict=None, labelpad=None, **csfont)

# trial
np.random.seed()
N=1000
x = np.random.uniform(-5.0, 5.0, N)
nbins = 100
plt.hist(x, nbins, normed=True, color='green',label="Rectangle") # Rectangle distribution

# font for legend
font = font_manager.FontProperties(family='Times New Roman',
                                   weight='bold',
                                   style='normal', size=10)
plt.legend(loc='upper left', prop=font) # legend

# plot options
plt.xticks([-5.0,-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0], **csfont)
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], **csfont)

# graph save & display
plt.savefig("histgram.png") # 1. file saving (1. should be before 2.)
plt.show() # 2. file showing (2. should be after 1.)

# paramters output
f.write(str(param[0])) # average
f1.write(str(param[1])) # standard deviation

# file close
f.close()
f1.close()

end()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
