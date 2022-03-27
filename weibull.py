# Weibull
# Coded by Takuro TOKUNAGA
# Last modified: September 11 2018

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
from scipy.stats import dweibull
from scipy.stats import beta
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

# graph
csfont = {'fontname':'Times New Roman'} # define font
plt.figure

temp1 = 0.45 # a
temp2 = 3.6 # c
temp3 = 0 # loc
temp4 = 3.5 # scale

x = np.linspace(-6,6, 100) # range of x axis for fitting curve
## self-defined distribution
#sns.kdeplot(data_1d_zsub, color = 'green', bw=0.75,label="bw: 0.75")
pdf_weib = exponweib.pdf(x, a=temp1, c=temp2, loc=temp3,scale=temp4) # gaussian fitting with parameters
plt.plot(x, pdf_weib, 'g-', label="Weibull (self-defined)") # display gaussian fitting with parameters

# random number generation
rand_weib = exponweib.rvs(a=temp1, c=temp2, loc=temp3,scale=temp4,size=1000)
#rand_weib = rand_weib-5
plt.hist(rand_weib, normed=True, histtype='stepfilled',alpha=.3, color='green', label='Weibull random') # not use density but use normed

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
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], **csfont)

# graph save & display
plt.savefig("histgram.png") # 1. file saving (1. should be before 2.)
plt.show()                  # 2. file showing (2. should be after 1.)

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
