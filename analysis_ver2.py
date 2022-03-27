# Output files analysis code
# Last modified: September 04 2018
# https://qiita.com/segavvy/items/51a515c19bcd29b13b7f
# Coded by Takuro TOKUNAGA

import math
import numpy as np
import cmath
import time
import sys
import pandas as pd
import glob
start = time.time()

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

txt_files = glob.glob("../surface/conductance/*.txt")
list = []
for f in txt_files:
    list.append(pd.read_csv(f))
df = pd.concat(list, axis=1)
df.to_csv("../surface/total_conductance.txt")

txt_files = glob.glob("../surface/area/*.txt")
list = []
for f in txt_files:
    list.append(pd.read_csv(f))
df = pd.concat(list, axis=1)
df.to_csv("../surface/total_area.txt")

txt_files = glob.glob("../surface/ratio/*.txt")
list = []
for f in txt_files:
    list.append(pd.read_csv(f))
df = pd.concat(list, axis=1)
df.to_csv("../surface/total_ratio.txt")

txt_files = glob.glob("../surface/cond_ac_ratio/*.txt")
list = []
for f in txt_files:
    list.append(pd.read_csv(f))
df = pd.concat(list, axis=1)
df.to_csv("../surface/total_cond_ac_ratio.txt")

txt_files = glob.glob("../surface/cond_nc_ratio/*.txt")
list = []
for f in txt_files:
    list.append(pd.read_csv(f))
df = pd.concat(list, axis=1)
df.to_csv("../surface/total_cond_nc_ratio.txt")

# dfdz (not used)
#txt_files = glob.glob("../surface/dfdz/*.txt")
#list = []
#for f in txt_files:
#    list.append(pd.read_csv(f))
#df = pd.concat(list, axis=1)
#df.to_csv("../surface/total_dfdz.txt")

## FE ##
txt_files = glob.glob("../surface/conductance_FE/*.txt")
list = []
for f in txt_files:
    list.append(pd.read_csv(f))
df = pd.concat(list, axis=1)
df.to_csv("../surface/total_conductance_FE.txt")

# calculation

# file close

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
