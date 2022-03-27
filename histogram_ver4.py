# Average & Std Calculator for surface roughness generator (sr_rand_gauss_ave.py)
# Surface roughness for Pt heater
# Last modified: June 03 2019
# Coded by Takuro TOKUNAGA

import numpy as np
import pandas as pd
# for graph
import matplotlib               # 1. matplotlib for mac
matplotlib.use('TkAgg')         # 2. matplotlib for mac
import matplotlib.pyplot as plt # 3. matplotlib for mac
import matplotlib.font_manager as font_manager
# for statistics
from scipy.stats import norm
import statistics
import pyper
r = pyper.R(use_pandas = "True")
import time
start = time.time()

# function: begin
def begin():
    print ("begin")

# function: end
def end():
    print ("end")

# unit conversion
ucnano = 1.0*np.power(10.,-9)

voltage1 = (-28.00) # [nm]
voltage2 = (-30.00) # [nm]
vfactor = (voltage1/voltage2)

# main
begin()

# file open
f1 = open('ave.txt', 'w') # write mode (outputs)
f2 = open('std.txt', 'w') # write mode (outputs)

## read raw data table, AFM_8107
data_2d = pd.read_csv("afm_snom_8107_21.txt", sep="\t", header=0, index_col=0)
row, col = data_2d.shape # row & column of matorix
data_1d_raw = np.ravel(data_2d) # 2d to 1d

# mean value of raw data
mean_value = statistics.mean(data_1d_raw) # mean
zshift = mean_value
data_1d_converted=np.zeros((row*col), dtype='float64') # sub: z shift subtract

# data manipulation: 0 point shift & unit conversion [nm]
for i in range(0,row*col):
    data_1d_converted[i] = vfactor*(data_1d_raw[i] - zshift)/ucnano # convert AFM raw data to Height data

max_converted = max(data_1d_converted) # max, maybe [nm]
min_converted = min(data_1d_converted) # min, maybe [nm]
mean_converted = statistics.mean(data_1d_converted) # mean, maybe [nm]
stdev_converted = statistics.stdev(data_1d_converted) # stdev
#print(str(max_value_zsub))
#print(str(min_value_zsub))
#print(str(mean_value_zsub))
#print(str(stdev_value_zsub))

## fitting parameters
param = norm.fit(data_1d_converted) # fitting paramter (average, standarad deviation)
# paramters output
f1.write(str(param[0])) # average
f2.write(str(param[1])) # standard deviation
# display
print("converted data average:{:.3f}".format(param[0]) + "[nm]")
print("converted data stdev:{:.3f}".format(param[1]) + "[nm]")

### for numerical bumpy start ###
# read numerical bumpy data
filenumber = 33
bumpy_2d = pd.read_csv("../surface/bumpy/data0.txt", sep=" ", header=None)
bumpy_2d.columns = ["x", "y", "h"]
row_numerical, col_numerical = bumpy_2d.shape # row & column of matorix
bumpy_1d_numerical = np.zeros(row_numerical*filenumber, dtype='float64')

# input data into tables
counter = 0
for i in range(0, filenumber):
    bumpy_2d = pd.read_csv("../surface/bumpy/data"+str(i)+".txt", sep=" ", header=None)
    bumpy_2d.columns = ["x", "y", "h"]
    for j in range(0, row_numerical):
        bumpy_1d_numerical[counter] = bumpy_2d.iat[j,2] # x line

        counter = counter + 1

# statistical values
min_numerical = min(bumpy_1d_numerical) # min [nm]
max_numerical = max(bumpy_1d_numerical) # max [nm]
mean_numerical = statistics.mean(bumpy_1d_numerical) # mean [nm]
stdev_numerical = statistics.stdev(bumpy_1d_numerical) # stdev [nm]

## fitting parameters
param_numerical = norm.fit(bumpy_1d_numerical) # fitting paramter (average, standarad deviation)

### for generated bumpy end ###

# graph display
csfont = {'fontname':'Times New Roman'} # define font
plt.figure
x = np.linspace(-6.0,6.0,100) # range of x axis for fitting curve

pdf_afm = norm.pdf(x,loc=param[0], scale=param[1]) # gaussian fitting with parameters
pdf_bumpy = norm.pdf(x,loc=param_numerical[0], scale=param_numerical[1]) # gaussian fitting with parameters

## for afm converted data histogram
weights = np.ones_like(data_1d_converted)/len(data_1d_converted)
bin_converted = int(np.sqrt(row*col)-8)

## for numerical bumpy histogram
weights_numerical = np.ones_like(bumpy_1d_numerical)/len(bumpy_1d_numerical)
bin_numerical = int(np.sqrt(row_numerical)-18)

plt.plot(x, pdf_afm, 'black', label="Gaussian distribution") # display gaussian fitting with parameters
plt.hist(data_1d_converted, normed=True, color = 'red', label = 'AFM surface roughness') # histgram of AFM data, default setting
plt.hist(bumpy_1d_numerical, weights=weights_numerical,bins=bin_numerical, alpha=0.3, color = 'blue', label = 'Numerical roughness') # histgram of AFM data, default setting

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
#stats.probplot(data_1d_converted, dist='norm', plot=ax) # AFM
#stats.probplot(bumpy_1d, dist='norm', plot=ax) # numerical
#plt.savefig("ppline.png") # 1. file saving (1. should be before 2.)
#plt.show()

#print(stats.shapiro(data_1d_converted))
#print(stats.shapiro(bumpy_1d))

### for origin graph ###
f5 = open('origin_histogram.txt', 'w') # read mode
f6 = open('origin_norm.txt', 'w') # read mode

## output, [0]: frequency, [1]: class value
ret = plt.hist(data_1d_converted, normed=True, color = 'red', label = 'AFM surface roughness') # histgram of AFM data, default setting
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
f5.close()
f6.close()
end()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
