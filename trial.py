# https://qiita.com/ynakayama/items/40272961bd43f63e6054
# https://qiita.com/yubais/items/bf9ce0a8fefdcc0b0c97

import numpy as np
import pandas as pd
import matplotlib               # 1. matplotlib for mac
matplotlib.use('TkAgg')         # 2. matplotlib for mac
import matplotlib.pyplot as plt # 3. matplotlib for mac
from scipy.stats import norm
from numpy.random import *

# data table
sample = norm.rvs(loc=100,scale=1,size=500)
# read by pandas

# fitting paramter (average, standarad deviation)
param = norm.fit(sample)
#print(param)

# graph
x = np.linspace(95,105,100) # range of x axis for fitting curve
pdf_fitted = norm.pdf(x,loc=param[0], scale=param[1]) # read fitting paramter
pdf = norm.pdf(x) # fitting curve
plt.figure
plt.title('Normal distribution') # graph title
plt.plot(x, pdf_fitted, 'r-', x,pdf, 'b-') # display fitting curve
plt.hist(sample, normed=1, alpha=.3) # histgram of original data
plt.savefig("image.png") # 1. file saving (1. should be before 2.)
plt.show() # 2. file showing (2. should be after 1.)

# random number generation
rnum = normal(param[0], param[1]) # random number
