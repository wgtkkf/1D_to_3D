# Average & Std Calculator for surface roughness generator (sr_rand_gauss_ave_tip.py)
# Surface roughness for Si tip
# Last modified: June 17 2019
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
import pylab as py
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
voltage1 = (-28.00) # [nm/V]
voltage2 = (-19.57) # [nm/V]
vfactor = (voltage1/voltage2)
vfactor = 1.0

# main
begin()

# file open
f1 = open('ave_tip.txt', 'w') # write mode (outputs)
f2 = open('std_tip.txt', 'w') # write mode (outputs)
f3 = open('converted afm data.txt', 'w') # write mode (outputs)

## read raw data table, afm_snom_####
#data_2d_raw = pd.read_csv("AFM_SNOM_8845_FWD.txt", sep="\t", header=0, index_col=0)
data_2d_raw = pd.read_csv("AFM_SNOM_8858.txt", sep="\t", header=0, index_col=0)

row, col = data_2d_raw.shape # row & column of matorix
data_1d_raw = np.ravel(data_2d_raw) # 2d to 1d

# mean value of raw data
mean_value = statistics.mean(data_1d_raw) # mean
zshift = mean_value
data_1d_converted=np.zeros((row*col), dtype='float64') # sub: z shift subtract

# raw data conversion: mean value shift & unit conversion [nm]
for i in range(0,row*col):
    data_1d_converted[i] = vfactor*(data_1d_raw[i] - zshift)/ucnano # [nm], convert AFM raw data to Height data
    # output converted data
    f3.write(str(i)) # [-]
    f3.write(' ')
    f3.write(str(data_1d_converted[i])) # [nm]
    f3.write('\n')

# quantities for converted data
max_converted = max(data_1d_converted) # max, maybe [nm]
min_converted = min(data_1d_converted) # min, maybe [nm]
mean_converted = statistics.mean(data_1d_converted) # mean, maybe [nm]
stdev_converted = statistics.stdev(data_1d_converted) # stdev
#print(str(max_converted))
#print(str(min_converted))
#print(str(mean_converted))
#print(str(stdev_converted))

## fitting parameters
param = norm.fit(data_1d_converted) # fitting paramter (average, standarad deviation)
# paramters output
f1.write(str(param[0])) # average [nm]
f2.write(str(param[1])) # standard deviation [nm]
# display
print("converted data average:{:.3f}".format(param[0]) + "[nm]")
print("converted data stdev:{:.3f}".format(param[1]) + "[nm]")
#print(str(param[0]))
#print(str(param[1]))

### for numerical bumpy: start ###
# read numerical bumpy data
filenumber = 33
bumpy_2d = pd.read_csv("../surface/bumpy_tip/data0.txt", sep=" ", header=None)
bumpy_2d.columns = ["x", "y", "h"]
row_numerical, col_numerical = bumpy_2d.shape # row & column of matorix
bumpy_1d_numerical = np.zeros(row_numerical*filenumber, dtype='float64')

# input data into tables
counter = 0
for i in range(0, filenumber):
    bumpy_2d = pd.read_csv("../surface/bumpy_tip/data"+str(i)+".txt", sep=" ", header=None)
    bumpy_2d.columns = ["x", "y", "h"]
    for j in range(0, row_numerical):
        bumpy_1d_numerical[counter] = bumpy_2d.iat[j,2] # x line

        counter = counter + 1

# statistical values
min_numerical = min(bumpy_1d_numerical) # min [nm]
max_numerical = max(bumpy_1d_numerical) # max [nm]
mean_numerical = statistics.mean(bumpy_1d_numerical) # mean [nm]
stdev_numerical = statistics.stdev(bumpy_1d_numerical) # stdev [nm]

# gaussian fitting parameters for numerical bumpy
param_numerical = norm.fit(bumpy_1d_numerical) # fitting paramter (average, standarad deviation)
#print(str(param_numerical[0]))
#print(str(param_numerical[1]))
### for numerical bumpy: end ###

# graph display
csfont = {'fontname':'Times New Roman'} # define font
plt.figure
x = np.linspace(-1.5,1.5,100) # range of x axis for fitting curve

# gaussian distribution for afm converted data
#pdf_afm = norm.pdf(x,loc=param[0], scale=param[1]) # gaussian fitting for afm converted data
sf = np.sqrt(2*np.pi*np.power(param[1],2.0)) # scale factor
#sf = np.sqrt(2*np.pi*np.power(param[1],2.0)) # scale factor
pdf_afm = (1-sf)*norm.pdf(x,loc=param[0], scale=param[1]) # gaussian fitting for afm converted data
# output weight


## for afm converted data histogram
weights = np.ones_like(data_1d_converted)/len(data_1d_converted)
bin_converted = int(np.sqrt(row*col)-10)

## for numerical bumpy histogram
weights_numerical = np.ones_like(bumpy_1d_numerical)/len(bumpy_1d_numerical)
bin_numerical = int(np.sqrt(row_numerical)-22)

# plot
plt.plot(x, pdf_afm, color = 'black', label="Gaussian distribution") # display gaussian fitting with parameters
plt.hist(data_1d_converted, weights=weights,bins=bin_converted, color = 'red', label = 'Measured roughness') # histgram of AFM data, default setting
plt.hist(bumpy_1d_numerical, weights=weights_numerical,bins=bin_numerical, alpha=0.3, color = 'blue', label = 'Numerical roughness') # histgram of AFM data, default setting

## graph information
#plt.title('Distribution', **csfont) # graph title
plt.xlabel('Tip Surface Roughness, H [nm]', fontdict=None, labelpad=None, **csfont)
plt.ylabel('Normalized Density [-]', fontdict=None, labelpad=None, **csfont)
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
## [0]: frequency, [1]: class value
### for drawing origin graph ###
f4 = open('origin_histogram_tip.txt', 'w') # write mode
f5 = open('origin_norm_tip.txt', 'w') # write mode

#ret = plt.hist(data_1d_converted, normed=True, color = 'red', label = 'Measured roughness') # histgram of AFM data, default setting
ret = plt.hist(data_1d_converted, weights=weights,bins=bin_converted, color = 'red', label = 'Measured roughness') # histgram of AFM data, default setting

frequency = ret[0].flatten()
for i in range(0,ret[0].shape[0]):
    f4.write(str(frequency[i])) # frequency
    f4.write('\n')

#f4.write(str(ret[0])) # frequency

f4.write('\n')
class_value = ret[1].flatten()
for i in range(0,ret[1].shape[0]):
    f4.write(str(class_value[i])) # class
    f4.write('\n')

#f4.write(str(ret[1])) # class value

## output for origin
## output of norm
normx = x.flatten()
normy = pdf_afm.flatten()
for i in range(0,x.shape[0]):
    f5.write(str(normx[i])) # class
    f5.write(' ')
    f5.write(str(normy[i])) # class
    f5.write('\n')

#f5.write(str(x)) # x
#f5.write('\n')
#f5.write('\n')
#f5.write(str(pdf_afm)) # norm

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
